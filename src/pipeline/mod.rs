pub mod worker;
pub mod writer;

use crate::error::{AppError, Result};
use crate::fastq::{InputReader, ReadOrPair};
use crate::output::summary_path;
use crate::stats::RunStats;
use crate::util::format_count;
use crossbeam_channel::{bounded, Receiver, Sender};
use std::path::Path;
use std::thread;
use std::time::{Duration, Instant};
use worker::{process_chunk, InputChunk, ProcessConfig, ProcessedChunk};
use writer::OrderedWriter;

/// Runtime options for writing outputs and threading.
pub struct PipelineOptions<'a> {
    pub input_files: &'a [&'a Path],
    pub out_dir: &'a Path,
    pub prefix: &'a str,
    pub gzip: bool,
    pub compression_level: u32,
    pub force: bool,
    pub threads: usize,
    pub chunk_reads: usize,
    pub summary: Option<&'a Path>,
    pub quiet: bool,
    /// 0 = no limit.
    pub max_reads: u64,
    /// Skip this many reads/pairs before counting toward max_reads.
    pub skip_reads: u64,
}

/// Run the full demultiplex pipeline.
pub fn run_pipeline(
    mut reader: InputReader,
    cfg: ProcessConfig,
    opts: PipelineOptions<'_>,
) -> Result<RunStats> {
    crate::util::validate_prefix(opts.prefix)?;
    let is_paired = match reader {
        InputReader::Paired(_) => true,
        InputReader::Single(_) => false,
    };
    let plan = crate::output::OutputPlan::build(&crate::output::OutputPlanParams {
        out_dir: opts.out_dir,
        prefix: opts.prefix,
        barcodes: &cfg.barcodes,
        is_paired,
        gzip: opts.gzip,
        discard_unassigned: cfg.discard_unassigned,
        counts_only: cfg.counts_only,
        summary: opts.summary,
    });
    crate::output::preflight_check(&plan, opts.input_files, opts.force)?;

    std::fs::create_dir_all(opts.out_dir)?;

    let threads = opts.threads.max(1);
    let chunk_reads = opts.chunk_reads.max(1);

    if threads == 1 {
        return run_serial(&mut reader, &cfg, &opts, chunk_reads);
    }

    let (job_tx, job_rx): (Sender<InputChunk>, Receiver<InputChunk>) = bounded(threads * 2);
    let (res_tx, res_rx): (
        Sender<Result<ProcessedChunk>>,
        Receiver<Result<ProcessedChunk>>,
    ) = bounded(threads * 2);

    let worker_cfg = cfg.clone();
    let mut handles = Vec::new();
    for _ in 0..threads {
        let rx = job_rx.clone();
        let tx = res_tx.clone();
        let wcfg = worker_cfg.clone();
        handles.push(thread::spawn(move || {
            while let Ok(chunk) = rx.recv() {
                let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    process_chunk(chunk, &wcfg)
                }));
                match result {
                    Ok(processed) => {
                        if tx.send(Ok(processed)).is_err() {
                            break;
                        }
                    }
                    Err(_) => {
                        let _ = tx.send(Err(AppError::WorkerFailure("worker panicked".into())));
                        break;
                    }
                }
            }
        }));
    }
    drop(res_tx);
    drop(job_rx);

    let max_reads = opts.max_reads;
    let skip_reads = opts.skip_reads;
    let reader_handle = thread::spawn(move || -> Result<u64> {
        skip_items(&mut reader, skip_reads)?;
        let mut chunk_id = 0u64;
        let mut n_read = 0u64;
        let mut batch = Vec::with_capacity(chunk_reads);
        while let Some(chunk) = next_chunk(
            &mut reader,
            &mut chunk_id,
            &mut batch,
            chunk_reads,
            max_reads,
            &mut n_read,
        )? {
            if job_tx.send(chunk).is_err() {
                return Err(AppError::WorkerFailure("workers disconnected".into()));
            }
        }
        Ok(chunk_id)
    });

    let mut writer = OrderedWriter::new(
        opts.out_dir,
        opts.prefix,
        opts.gzip,
        opts.compression_level,
        opts.force,
    );
    let mut received = 0u64;
    let start = Instant::now();
    let mut last_report = Instant::now();

    loop {
        match res_rx.recv() {
            Ok(Ok(chunk)) => {
                received += 1;
                writer.push(chunk)?;
                if should_report(opts.quiet, &mut last_report) {
                    report_progress(&writer.stats, start);
                }
            }
            Ok(Err(e)) => return Err(e),
            Err(_) => break,
        }
    }

    let n_chunks = reader_handle
        .join()
        .map_err(|_| AppError::WorkerFailure("reader thread panicked".into()))??;

    for h in handles {
        let _ = h.join();
    }

    if received != n_chunks {
        return Err(AppError::WorkerFailure(format!(
            "received {received} processed chunks, expected {n_chunks} (worker failure or dropped chunks)"
        )));
    }

    let stats = writer.finish()?;
    write_summary(&stats, opts.out_dir, opts.prefix, opts.summary, opts.quiet)?;
    Ok(stats)
}

fn run_serial(
    reader: &mut InputReader,
    cfg: &ProcessConfig,
    opts: &PipelineOptions<'_>,
    chunk_reads: usize,
) -> Result<RunStats> {
    let mut writer = OrderedWriter::new(
        opts.out_dir,
        opts.prefix,
        opts.gzip,
        opts.compression_level,
        opts.force,
    );
    let mut chunk_id = 0u64;
    let mut n_read = 0u64;
    let mut batch = Vec::with_capacity(chunk_reads);
    let start = Instant::now();
    let mut last_report = Instant::now();
    skip_items(reader, opts.skip_reads)?;

    while let Some(chunk) = next_chunk(
        reader,
        &mut chunk_id,
        &mut batch,
        chunk_reads,
        opts.max_reads,
        &mut n_read,
    )? {
        let processed = process_chunk(chunk, cfg);
        writer.push(processed)?;
        if should_report(opts.quiet, &mut last_report) {
            report_progress(&writer.stats, start);
        }
    }

    let stats = writer.finish()?;
    write_summary(&stats, opts.out_dir, opts.prefix, opts.summary, opts.quiet)?;
    Ok(stats)
}

fn write_summary(
    stats: &RunStats,
    out_dir: &Path,
    prefix: &str,
    summary: Option<&Path>,
    quiet: bool,
) -> Result<()> {
    let summary_out = summary
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| summary_path(out_dir, prefix));
    stats.write_tsv(&summary_out)?;
    if !quiet {
        stats.print_summary();
        eprintln!("Summary written to {}", summary_out.display());
    }
    Ok(())
}

fn next_chunk(
    reader: &mut InputReader,
    chunk_id: &mut u64,
    batch: &mut Vec<ReadOrPair>,
    chunk_reads: usize,
    max_reads: u64,
    n_read: &mut u64,
) -> Result<Option<InputChunk>> {
    loop {
        if max_reads > 0 && *n_read >= max_reads {
            break;
        }
        match reader.next_item()? {
            Some(item) => {
                *n_read += 1;
                batch.push(item);
                if batch.len() >= chunk_reads {
                    break;
                }
            }
            None => break,
        }
    }
    if batch.is_empty() {
        return Ok(None);
    }
    let chunk = InputChunk {
        id: *chunk_id,
        records: std::mem::replace(batch, Vec::with_capacity(chunk_reads)),
    };
    *chunk_id += 1;
    Ok(Some(chunk))
}

fn skip_items(reader: &mut InputReader, skip: u64) -> Result<()> {
    for _ in 0..skip {
        if reader.next_item()?.is_none() {
            break;
        }
    }
    Ok(())
}

fn should_report(quiet: bool, last_report: &mut Instant) -> bool {
    if quiet {
        return false;
    }
    if last_report.elapsed() >= Duration::from_secs(2) {
        *last_report = Instant::now();
        true
    } else {
        false
    }
}

fn report_progress(stats: &RunStats, start: Instant) {
    let elapsed = start.elapsed().as_secs_f64().max(1e-9);
    let n = stats.total_reads;
    if n == 0 {
        return;
    }
    let rate = n as f64 / elapsed;
    let assigned_pct = 100.0 * stats.assigned as f64 / n as f64;
    eprintln!(
        "processed {} pairs | {:.1} k pairs/s | assigned {} ({:.1}%) | {:.0}s",
        format_count(n),
        rate / 1000.0,
        format_count(stats.assigned),
        assigned_pct,
        elapsed
    );
}
