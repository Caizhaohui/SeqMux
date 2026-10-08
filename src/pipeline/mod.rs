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
    if opts.threads < 1 {
        return Err(AppError::Cli("threads must be at least 1".into()));
    }
    if opts.chunk_reads < 1 {
        return Err(AppError::Cli("chunk_reads must be at least 1".into()));
    }
    if !(1..=9).contains(&opts.compression_level) {
        return Err(AppError::Cli(
            "compression_level must be between 1 and 9".into(),
        ));
    }
    if !cfg.adapter_error_rate.is_finite() || !(0.0..=1.0).contains(&cfg.adapter_error_rate) {
        return Err(AppError::Cli(
            "adapter_error_rate must be a finite number between 0.0 and 1.0".into(),
        ));
    }
    if cfg.min_adapter_overlap < 1 {
        return Err(AppError::Cli(
            "min_adapter_overlap must be at least 1".into(),
        ));
    }

    crate::util::validate_prefix(opts.prefix)?;
    let is_paired = match reader {
        InputReader::Paired(_) | InputReader::ConcurrentPaired(_) => true,
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

    let threads = opts.threads;
    let chunk_reads = opts.chunk_reads;

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
        match reader {
            InputReader::ConcurrentPaired(mut cr) => {
                let mut n_chunks = 0u64;
                while let Some((chunk_id, pairs)) = cr.next_paired_chunk()? {
                    n_chunks += 1;
                    if job_tx
                        .send(InputChunk {
                            id: chunk_id,
                            records: pairs,
                        })
                        .is_err()
                    {
                        return Err(AppError::WorkerFailure("workers disconnected".into()));
                    }
                }
                cr.finish()?;
                Ok(n_chunks)
            }
            mut other_reader => {
                skip_items(&mut other_reader, skip_reads)?;
                let mut chunk_id = 0u64;
                let mut n_read = 0u64;
                let mut batch = Vec::with_capacity(chunk_reads);
                while let Some(chunk) = next_chunk(
                    &mut other_reader,
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
            }
        }
    });

    let sample_labels: Vec<String> = cfg
        .barcodes
        .samples
        .iter()
        .map(|s| s.sanitized_label.clone())
        .collect();
    let mut writer = OrderedWriter::new(
        opts.out_dir,
        opts.prefix,
        opts.gzip,
        opts.compression_level,
        opts.force,
        sample_labels,
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
    let sample_labels: Vec<String> = cfg
        .barcodes
        .samples
        .iter()
        .map(|s| s.sanitized_label.clone())
        .collect();
    let mut writer = OrderedWriter::new(
        opts.out_dir,
        opts.prefix,
        opts.gzip,
        opts.compression_level,
        opts.force,
        sample_labels,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::barcode::{BarcodeConfig, DemuxMode, OrientationMode};
    use std::path::Path;
    use tempfile::NamedTempFile;

    fn dummy_cfg() -> ProcessConfig {
        ProcessConfig {
            barcodes: BarcodeConfig {
                samples: Vec::new(),
                mode: DemuxMode::SingleBarcode,
                mismatches_1: 0,
                mismatches_2: 0,
                min_mismatch_delta: 0,
                fast_exact_8bp_pe: None,
            },
            keep_barcodes: false,
            discard_unassigned: false,
            min_length: 0,
            quality_cutoff_5: 0,
            quality_cutoff_3: 0,
            nextseq: false,
            adapter_r1: None,
            adapter_r2: None,
            adapter_error_rate: 0.1,
            min_adapter_overlap: 3,
            tso: None,
            phred_offset: 33,
            orientation: OrientationMode::Both,
            canonicalize: true,
            counts_only: false,
        }
    }

    fn dummy_opts<'a>(out_dir: &'a Path, input_files: &'a [&'a Path]) -> PipelineOptions<'a> {
        PipelineOptions {
            input_files,
            out_dir,
            prefix: "test",
            gzip: false,
            compression_level: 6,
            force: true,
            threads: 1,
            chunk_reads: 100,
            summary: None,
            quiet: true,
            max_reads: 0,
            skip_reads: 0,
        }
    }

    fn write_fq(path: &Path) {
        std::fs::write(path, "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    }

    #[test]
    fn test_run_pipeline_rejects_zero_threads() {
        let f = NamedTempFile::new().unwrap();
        let input_path = f.path();
        write_fq(input_path);
        let input_files = [input_path];
        let out_dir = tempfile::tempdir().unwrap();
        let reader = InputReader::single(input_path).unwrap();
        let cfg = dummy_cfg();
        let mut opts = dummy_opts(out_dir.path(), &input_files);
        opts.threads = 0;

        match run_pipeline(reader, cfg, opts) {
            Err(e) => assert!(e.to_string().contains("threads must be at least 1")),
            Ok(_) => panic!("expected error for threads=0"),
        }
    }

    #[test]
    fn test_run_pipeline_rejects_zero_chunk_reads() {
        let f = NamedTempFile::new().unwrap();
        let input_path = f.path();
        write_fq(input_path);
        let input_files = [input_path];
        let out_dir = tempfile::tempdir().unwrap();
        let reader = InputReader::single(input_path).unwrap();
        let cfg = dummy_cfg();
        let mut opts = dummy_opts(out_dir.path(), &input_files);
        opts.chunk_reads = 0;

        match run_pipeline(reader, cfg, opts) {
            Err(e) => assert!(e.to_string().contains("chunk_reads must be at least 1")),
            Ok(_) => panic!("expected error for chunk_reads=0"),
        }
    }

    #[test]
    fn test_run_pipeline_rejects_invalid_compression_level() {
        let f = NamedTempFile::new().unwrap();
        let input_path = f.path();
        write_fq(input_path);
        let input_files = [input_path];
        let out_dir = tempfile::tempdir().unwrap();
        let cfg = dummy_cfg();

        let reader0 = InputReader::single(input_path).unwrap();
        let mut opts0 = dummy_opts(out_dir.path(), &input_files);
        opts0.compression_level = 0;
        match run_pipeline(reader0, cfg.clone(), opts0) {
            Err(e) => assert!(e
                .to_string()
                .contains("compression_level must be between 1 and 9")),
            Ok(_) => panic!("expected error for compression_level=0"),
        }

        let reader10 = InputReader::single(input_path).unwrap();
        let mut opts10 = dummy_opts(out_dir.path(), &input_files);
        opts10.compression_level = 10;
        match run_pipeline(reader10, cfg, opts10) {
            Err(e) => assert!(e
                .to_string()
                .contains("compression_level must be between 1 and 9")),
            Ok(_) => panic!("expected error for compression_level=10"),
        }
    }

    #[test]
    fn test_run_pipeline_rejects_invalid_adapter_settings() {
        let f = NamedTempFile::new().unwrap();
        let input_path = f.path();
        write_fq(input_path);
        let input_files = [input_path];
        let out_dir = tempfile::tempdir().unwrap();

        // Negative adapter_error_rate
        let mut cfg_neg = dummy_cfg();
        cfg_neg.adapter_error_rate = -0.1;
        let reader_neg = InputReader::single(input_path).unwrap();
        match run_pipeline(
            reader_neg,
            cfg_neg,
            dummy_opts(out_dir.path(), &input_files),
        ) {
            Err(e) => assert!(e.to_string().contains("adapter_error_rate")),
            Ok(_) => panic!("expected error for negative adapter_error_rate"),
        }

        // NaN adapter_error_rate
        let mut cfg_nan = dummy_cfg();
        cfg_nan.adapter_error_rate = f32::NAN;
        let reader_nan = InputReader::single(input_path).unwrap();
        match run_pipeline(
            reader_nan,
            cfg_nan,
            dummy_opts(out_dir.path(), &input_files),
        ) {
            Err(e) => assert!(e.to_string().contains("adapter_error_rate")),
            Ok(_) => panic!("expected error for NaN adapter_error_rate"),
        }

        // Zero min_adapter_overlap
        let mut cfg_zero_ov = dummy_cfg();
        cfg_zero_ov.min_adapter_overlap = 0;
        let reader_ov = InputReader::single(input_path).unwrap();
        match run_pipeline(
            reader_ov,
            cfg_zero_ov,
            dummy_opts(out_dir.path(), &input_files),
        ) {
            Err(e) => assert!(e.to_string().contains("min_adapter_overlap")),
            Ok(_) => panic!("expected error for zero min_adapter_overlap"),
        }
    }
}
