use crate::barcode::{load_barcodes_csv, TsoPattern};
use crate::error::{AppError, Result};
use crate::fastq::{validate_fastq, InputReader};
use crate::pipeline::run_pipeline;
use crate::pipeline::worker::ProcessConfig;
use crate::pipeline::PipelineOptions;
use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "seqmux",
    version,
    about = "SeqMux — fast, portable FASTQ demultiplexer",
    long_about = "SeqMux demultiplexes single-end and paired-end FASTQ files using a \
sample barcode table (SampleNumber, Barcode1, Barcode2). Supports quality/adapter \
trimming and multi-threaded streaming I/O.\n\n\
Barcode table format is SeqMux-native (header required). Ultraplex CSV is not supported."
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Demultiplex FASTQ reads by barcode
    Demux(Box<DemuxArgs>),
    /// Validate FASTQ structure and pair consistency
    Validate(ValidateArgs),
}

#[derive(Parser, Debug)]
pub struct DemuxArgs {
    /// Input FASTQ (.fastq/.fq, optionally .gz)
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Paired-end R2 FASTQ
    #[arg(short = 'I', long = "input2", visible_alias = "i2")]
    pub input2: Option<PathBuf>,

    /// Sample barcode table CSV (header: SampleNumber,Barcode1,Barcode2,…)
    #[arg(short = 'b', long = "barcodes")]
    pub barcodes: PathBuf,

    /// Output directory
    #[arg(short = 'o', long = "out-dir", default_value = ".")]
    pub out_dir: PathBuf,

    /// Output filename prefix
    #[arg(short = 'p', long = "prefix", default_value = "seqmux")]
    pub prefix: String,

    /// Gzip compression level (1-9)
    #[arg(long = "compression-level", default_value_t = 6)]
    pub compression_level: u32,

    /// Write uncompressed .fastq instead of .fastq.gz
    #[arg(long = "no-gzip")]
    pub no_gzip: bool,

    /// Discard unassigned reads (do not write unassigned file)
    #[arg(long = "discard-unassigned")]
    pub discard_unassigned: bool,

    /// Allowed mismatches for Barcode1
    #[arg(
        long = "mismatches-1",
        visible_alias = "m1",
        alias = "mismatches-5",
        default_value_t = 0
    )]
    pub mismatches_1: usize,

    /// Allowed mismatches for Barcode2
    #[arg(
        long = "mismatches-2",
        visible_alias = "m2",
        alias = "mismatches-3",
        default_value_t = 0
    )]
    pub mismatches_2: usize,

    /// Keep barcode bases in the sequence (do not trim)
    #[arg(long = "keep-barcodes", visible_alias = "kbc")]
    pub keep_barcodes: bool,

    /// TSO pattern (N=UMI, I=ignore/trim)
    #[arg(long = "tso-pattern")]
    pub tso_pattern: Option<String>,

    /// 3' adapter sequence for R1
    #[arg(short = 'a', long = "adapter-r1")]
    pub adapter_r1: Option<String>,

    /// 3' adapter sequence for R2
    #[arg(long = "adapter-r2")]
    pub adapter_r2: Option<String>,

    /// 3' quality cutoff (Phred)
    #[arg(short = 'q', long = "quality-cutoff-3", default_value_t = 0)]
    pub quality_cutoff_3: u8,

    /// 5' quality cutoff (Phred)
    #[arg(long = "quality-cutoff-5", visible_alias = "q5", default_value_t = 0)]
    pub quality_cutoff_5: u8,

    /// NextSeq-style quality trimming
    #[arg(long = "nextseq")]
    pub nextseq: bool,

    /// Minimum read length after trimming
    #[arg(short = 'l', long = "min-length", default_value_t = 0)]
    pub min_length: usize,

    /// Minimum adapter overlap for trimming
    #[arg(long = "min-adapter-overlap", default_value_t = 3)]
    pub min_adapter_overlap: usize,

    /// Max error rate for adapter matching
    #[arg(long = "adapter-error-rate", default_value_t = 0.1)]
    pub adapter_error_rate: f32,

    /// Worker threads
    #[arg(short = 't', long = "threads", default_value_t = 4)]
    pub threads: usize,

    /// Reads per processing chunk
    #[arg(long = "chunk-reads", default_value_t = 4096)]
    pub chunk_reads: usize,

    /// Summary TSV path (default: <out-dir>/<prefix>.summary.tsv)
    #[arg(long = "summary")]
    pub summary: Option<PathBuf>,

    /// Overwrite existing output files
    #[arg(long = "force")]
    pub force: bool,

    /// Suppress progress and summary on stderr
    #[arg(long = "quiet")]
    pub quiet: bool,

    /// Log level
    #[arg(long = "log-level", default_value = "info")]
    pub log_level: LogLevel,
}

#[derive(Parser, Debug)]
pub struct ValidateArgs {
    /// Input FASTQ
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Optional R2 FASTQ
    #[arg(short = 'I', long = "input2")]
    pub input2: Option<PathBuf>,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum LogLevel {
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

impl LogLevel {
    pub fn as_str(self) -> &'static str {
        match self {
            LogLevel::Error => "error",
            LogLevel::Warn => "warn",
            LogLevel::Info => "info",
            LogLevel::Debug => "debug",
            LogLevel::Trace => "trace",
        }
    }
}

pub fn run(cli: Cli) -> Result<()> {
    match cli.command {
        Commands::Demux(args) => run_demux(*args),
        Commands::Validate(args) => run_validate(args),
    }
}

fn run_demux(args: DemuxArgs) -> Result<()> {
    if !args.quiet {
        env_logger::Builder::from_env(
            env_logger::Env::default().default_filter_or(args.log_level.as_str()),
        )
        .format_timestamp_secs()
        .init();
    }

    if !(1..=9).contains(&args.compression_level) {
        return Err(AppError::Cli(
            "--compression-level must be between 1 and 9".into(),
        ));
    }

    let barcodes = load_barcodes_csv(&args.barcodes, args.mismatches_1, args.mismatches_2)?;

    let tso = match args.tso_pattern {
        Some(ref s) => Some(TsoPattern::parse(s)?),
        None => None,
    };

    let cfg = ProcessConfig {
        barcodes,
        keep_barcodes: args.keep_barcodes,
        discard_unassigned: args.discard_unassigned,
        min_length: args.min_length,
        quality_cutoff_5: args.quality_cutoff_5,
        quality_cutoff_3: args.quality_cutoff_3,
        nextseq: args.nextseq,
        adapter_r1: args.adapter_r1.map(|s| s.into_bytes()),
        adapter_r2: args.adapter_r2.map(|s| s.into_bytes()),
        adapter_error_rate: args.adapter_error_rate,
        min_adapter_overlap: args.min_adapter_overlap,
        tso,
        phred_offset: 33,
    };

    let reader = if let Some(ref i2) = args.input2 {
        InputReader::paired(&args.input, i2)?
    } else {
        InputReader::single(&args.input)?
    };

    if !args.quiet {
        eprintln!("SeqMux demux");
        eprintln!("  input : {}", args.input.display());
        if let Some(ref i2) = args.input2 {
            eprintln!("  input2: {}", i2.display());
        }
        eprintln!("  barcodes: {}", args.barcodes.display());
        eprintln!("  samples : {}", cfg.barcodes.samples.len());
        eprintln!("  mode    : {:?}", cfg.barcodes.mode);
        eprintln!("  out-dir : {}", args.out_dir.display());
        eprintln!("  threads : {}", args.threads);
    }

    run_pipeline(
        reader,
        cfg,
        PipelineOptions {
            out_dir: &args.out_dir,
            prefix: &args.prefix,
            gzip: !args.no_gzip,
            compression_level: args.compression_level,
            force: args.force,
            threads: args.threads,
            chunk_reads: args.chunk_reads,
            summary: args.summary.as_deref(),
            quiet: args.quiet,
        },
    )?;

    Ok(())
}

fn run_validate(args: ValidateArgs) -> Result<()> {
    let stats = validate_fastq(&args.input, args.input2.as_ref())?;
    println!("records\t{}", stats.records);
    println!("bases\t{}", stats.bases);
    println!(
        "min_length\t{}",
        stats
            .min_length
            .map(|n| n.to_string())
            .unwrap_or_else(|| "NA".into())
    );
    println!(
        "max_length\t{}",
        stats
            .max_length
            .map(|n| n.to_string())
            .unwrap_or_else(|| "NA".into())
    );
    println!("paired\t{}", stats.paired);
    println!("status\tOK");
    Ok(())
}
