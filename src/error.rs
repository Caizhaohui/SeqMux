use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum AppError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("FASTQ format error: {0}")]
    FastqFormat(String),

    #[error("barcode config error: {0}")]
    BarcodeConfig(String),

    #[error("paired FASTQ records are out of sync\n  R1 record {index}: {r1}\n  R2 record {index}: {r2}")]
    InvalidPair { index: u64, r1: String, r2: String },

    #[error("invalid base in barcode: {0}")]
    InvalidBase(String),

    #[error("output conflict: {path} already exists (use --force to overwrite)")]
    OutputConflict { path: PathBuf },

    #[error("worker failure: {0}")]
    WorkerFailure(String),

    #[error("CLI error: {0}")]
    Cli(String),

    #[error("{0}")]
    Other(String),
}

pub type Result<T> = std::result::Result<T, AppError>;

impl From<anyhow::Error> for AppError {
    fn from(err: anyhow::Error) -> Self {
        AppError::Other(err.to_string())
    }
}

impl From<csv::Error> for AppError {
    fn from(err: csv::Error) -> Self {
        AppError::BarcodeConfig(err.to_string())
    }
}
