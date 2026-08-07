use crate::barcode::SampleKey;
use std::path::{Path, PathBuf};

/// Output file identity.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct OutputKey {
    pub sample: SampleKey,
    pub mate: Mate,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Mate {
    Single,
    R1,
    R2,
}

/// Build output path for a sample.
pub fn output_path(
    out_dir: &Path,
    prefix: &str,
    sample: &SampleKey,
    mate: Mate,
    gzip: bool,
) -> PathBuf {
    let label = sample.label();
    let ext = if gzip { "fastq.gz" } else { "fastq" };
    let name = match mate {
        Mate::Single => format!("{prefix}_{label}.{ext}"),
        Mate::R1 => format!("{prefix}_{label}_R1.{ext}"),
        Mate::R2 => format!("{prefix}_{label}_R2.{ext}"),
    };
    out_dir.join(name)
}

/// Summary TSV path.
pub fn summary_path(out_dir: &Path, prefix: &str) -> PathBuf {
    out_dir.join(format!("{prefix}.summary.tsv"))
}
