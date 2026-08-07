//! SeqMux — portable FASTQ demultiplexer.
//!
//! Core library modules for barcode demultiplexing, quality/adapter trimming,
//! and multi-threaded FASTQ I/O.

pub mod barcode;
pub mod cli;
pub mod error;
pub mod fastq;
pub mod output;
pub mod pipeline;
pub mod stats;
pub mod trim;
pub mod util;

pub use error::{AppError, Result};
