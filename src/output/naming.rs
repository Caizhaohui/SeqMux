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

use crate::barcode::BarcodeConfig;
use crate::error::{AppError, Result};
use crate::util::paths_point_to_same_file;
use std::collections::HashSet;

/// Parameters for constructing an OutputPlan.
#[derive(Debug, Clone)]
pub struct OutputPlanParams<'a> {
    pub out_dir: &'a Path,
    pub prefix: &'a str,
    pub barcodes: &'a BarcodeConfig,
    pub is_paired: bool,
    pub gzip: bool,
    pub discard_unassigned: bool,
    pub counts_only: bool,
    pub summary: Option<&'a Path>,
}

/// Planned output paths for a demultiplexing run.
#[derive(Debug, Clone)]
pub struct OutputPlan {
    pub fastq_files: Vec<PathBuf>,
    pub summary_file: PathBuf,
}

impl OutputPlan {
    pub fn build(params: &OutputPlanParams<'_>) -> Self {
        let summary_file = params
            .summary
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| summary_path(params.out_dir, params.prefix));

        if params.counts_only {
            return Self {
                fastq_files: Vec::new(),
                summary_file,
            };
        }

        let mut fastq_files = Vec::new();
        for sample in &params.barcodes.samples {
            let key = SampleKey::Named(sample.name.clone());
            if params.is_paired {
                fastq_files.push(output_path(
                    params.out_dir,
                    params.prefix,
                    &key,
                    Mate::R1,
                    params.gzip,
                ));
                fastq_files.push(output_path(
                    params.out_dir,
                    params.prefix,
                    &key,
                    Mate::R2,
                    params.gzip,
                ));
            } else {
                fastq_files.push(output_path(
                    params.out_dir,
                    params.prefix,
                    &key,
                    Mate::Single,
                    params.gzip,
                ));
            }
        }

        if !params.discard_unassigned {
            let key = SampleKey::Unassigned;
            if params.is_paired {
                fastq_files.push(output_path(
                    params.out_dir,
                    params.prefix,
                    &key,
                    Mate::R1,
                    params.gzip,
                ));
                fastq_files.push(output_path(
                    params.out_dir,
                    params.prefix,
                    &key,
                    Mate::R2,
                    params.gzip,
                ));
            } else {
                fastq_files.push(output_path(
                    params.out_dir,
                    params.prefix,
                    &key,
                    Mate::Single,
                    params.gzip,
                ));
            }
        }

        Self {
            fastq_files,
            summary_file,
        }
    }

    /// All files that SeqMux will write for this run.
    pub fn all_output_files(&self) -> Vec<&Path> {
        let mut list: Vec<&Path> = self.fastq_files.iter().map(|p| p.as_path()).collect();
        list.push(&self.summary_file);
        list
    }
}

/// Perform preflight validation of all planned output files.
pub fn preflight_check(plan: &OutputPlan, input_files: &[&Path], force: bool) -> Result<()> {
    // 1. Check collisions among planned outputs
    let mut seen_canonical = HashSet::new();
    for path in plan.all_output_files() {
        let norm = crate::util::normalize_path_for_compare(path);
        if !seen_canonical.insert(norm) {
            return Err(AppError::OutputConflict {
                path: path.to_path_buf(),
            });
        }
    }

    // 2. Check collisions with input files (fatal even with --force)
    for out_path in plan.all_output_files() {
        for &in_path in input_files {
            if paths_point_to_same_file(out_path, in_path) {
                return Err(AppError::Cli(format!(
                    "output file '{}' would overwrite input file '{}'",
                    out_path.display(),
                    in_path.display()
                )));
            }
        }
    }

    // 3. Check for existing outputs if not --force
    if !force {
        let mut existing = Vec::new();
        for path in plan.all_output_files() {
            if path.exists() {
                existing.push(path.display().to_string());
            }
        }
        if !existing.is_empty() {
            return Err(AppError::OutputConflict {
                path: PathBuf::from(existing.join(", ")),
            });
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::barcode::parse_barcodes_csv;
    use tempfile::tempdir;

    #[test]
    fn preflight_detects_existing_file_without_force() {
        let dir = tempdir().unwrap();
        let out_dir = dir.path().join("out");
        std::fs::create_dir_all(&out_dir).unwrap();
        let existing = out_dir.join("seqmux_A1_R1.fastq.gz");
        std::fs::write(&existing, b"test").unwrap();

        let csv = "SampleNumber,Barcode1,Barcode2\nA1,AAGTCCAA,GGAGTACT\n";
        let cfg = parse_barcodes_csv(csv, 0, 0).unwrap();
        let plan = OutputPlan::build(&OutputPlanParams {
            out_dir: &out_dir,
            prefix: "seqmux",
            barcodes: &cfg,
            is_paired: true,
            gzip: true,
            discard_unassigned: false,
            counts_only: false,
            summary: None,
        });

        let in_file = dir.path().join("in.fastq");
        assert!(preflight_check(&plan, &[&in_file], false).is_err());
        assert!(preflight_check(&plan, &[&in_file], true).is_ok());
    }

    #[test]
    fn preflight_detects_output_overwriting_input() {
        let dir = tempdir().unwrap();
        let out_dir = dir.path().join("out");
        std::fs::create_dir_all(&out_dir).unwrap();

        let csv = "SampleNumber,Barcode1,Barcode2\nA1,AAGTCCAA,GGAGTACT\n";
        let cfg = parse_barcodes_csv(csv, 0, 0).unwrap();
        let plan = OutputPlan::build(&OutputPlanParams {
            out_dir: &out_dir,
            prefix: "seqmux",
            barcodes: &cfg,
            is_paired: true,
            gzip: true,
            discard_unassigned: false,
            counts_only: false,
            summary: None,
        });

        let input_matching_output = out_dir.join("seqmux_A1_R1.fastq.gz");
        // Even with force = true, overwriting input is fatal
        let err = preflight_check(&plan, &[&input_matching_output], true).unwrap_err();
        assert!(err.to_string().contains("would overwrite input file"));
    }

    #[test]
    fn preflight_detects_summary_collision_with_fastq() {
        let dir = tempdir().unwrap();
        let out_dir = dir.path().join("out");

        let csv = "SampleNumber,Barcode1\nA1,AAGTCCAA\n";
        let cfg = parse_barcodes_csv(csv, 0, 0).unwrap();
        // Set summary path to collide with A1 output file
        let colliding_summary = out_dir.join("seqmux_A1.fastq");
        let plan = OutputPlan::build(&OutputPlanParams {
            out_dir: &out_dir,
            prefix: "seqmux",
            barcodes: &cfg,
            is_paired: false,
            gzip: false,
            discard_unassigned: false,
            counts_only: false,
            summary: Some(&colliding_summary),
        });

        let in_file = dir.path().join("in.fastq");
        let err = preflight_check(&plan, &[&in_file], true).unwrap_err();
        assert!(matches!(err, AppError::OutputConflict { .. }));
    }
}
