pub mod worker;
pub mod writer;

use crate::error::{AppError, Result};
use crate::fastq::InputReader;
use crate::output::summary_path;
use crate::stats::RunStats;
use crossbeam_channel::{bounded, Receiver, Sender};
use std::path::Path;
use std::thread;
use worker::{process_chunk, InputChunk, ProcessConfig, ProcessedChunk};
use writer::OrderedWriter;

/// Runtime options for writing outputs and threading.
pub struct PipelineOptions<'a> {
    pub out_dir: &'a Path,
    pub prefix: &'a str,
    pub gzip: bool,
    pub compression_level: u32,
    pub force: bool,
    pub threads: usize,
    pub chunk_reads: usize,
    pub summary: Option<&'a Path>,
    pub quiet: bool,
}

/// Run the full demultiplex pipeline.
pub fn run_pipeline(
    mut reader: InputReader,
    cfg: ProcessConfig,
    opts: PipelineOptions<'_>,
) -> Result<RunStats> {
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

    let reader_handle = thread::spawn(move || -> Result<u64> {
        let mut chunk_id = 0u64;
        let mut batch = Vec::with_capacity(chunk_reads);
        loop {
            match reader.next_item()? {
                Some(item) => {
                    batch.push(item);
                    if batch.len() >= chunk_reads {
                        let chunk = InputChunk {
                            id: chunk_id,
                            records: std::mem::replace(&mut batch, Vec::with_capacity(chunk_reads)),
                        };
                        chunk_id += 1;
                        if job_tx.send(chunk).is_err() {
                            return Err(AppError::WorkerFailure("workers disconnected".into()));
                        }
                    }
                }
                None => {
                    if !batch.is_empty() {
                        let chunk = InputChunk {
                            id: chunk_id,
                            records: batch,
                        };
                        chunk_id += 1;
                        if job_tx.send(chunk).is_err() {
                            return Err(AppError::WorkerFailure("workers disconnected".into()));
                        }
                    }
                    break;
                }
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

    loop {
        match res_rx.recv() {
            Ok(Ok(chunk)) => {
                received += 1;
                if !opts.quiet && received % 50 == 0 {
                    eprintln!(
                        "processed ~{} chunks | assigned {}",
                        received, writer.stats.assigned
                    );
                }
                writer.push(chunk)?;
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
        log::warn!(
            "received {received} processed chunks, expected {n_chunks} (check for worker errors)"
        );
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
    let mut batch = Vec::with_capacity(chunk_reads);

    loop {
        match reader.next_item()? {
            Some(item) => {
                batch.push(item);
                if batch.len() >= chunk_reads {
                    let chunk = InputChunk {
                        id: chunk_id,
                        records: std::mem::replace(&mut batch, Vec::with_capacity(chunk_reads)),
                    };
                    chunk_id += 1;
                    let processed = process_chunk(chunk, cfg);
                    writer.push(processed)?;
                }
            }
            None => {
                if !batch.is_empty() {
                    let chunk = InputChunk {
                        id: chunk_id,
                        records: batch,
                    };
                    let processed = process_chunk(chunk, cfg);
                    writer.push(processed)?;
                }
                break;
            }
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
