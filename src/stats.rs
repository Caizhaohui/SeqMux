use crate::barcode::SampleKey;
use crate::util::format_count;
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::io::Write;
use std::path::Path;

#[derive(Debug, Default, Clone)]
pub struct ChunkStats {
    pub total_reads: u64,
    pub quality_trimmed: u64,
    pub adapter_trimmed: u64,
    pub assigned: u64,
    pub unassigned: u64,
    pub ambiguous: u64,
    pub too_short: u64,
    pub five_prime_matched_three_prime_missing: u64,
    pub per_sample: BTreeMap<String, u64>,
}

impl ChunkStats {
    pub fn merge(&mut self, other: &ChunkStats) {
        self.total_reads += other.total_reads;
        self.quality_trimmed += other.quality_trimmed;
        self.adapter_trimmed += other.adapter_trimmed;
        self.assigned += other.assigned;
        self.unassigned += other.unassigned;
        self.ambiguous += other.ambiguous;
        self.too_short += other.too_short;
        self.five_prime_matched_three_prime_missing += other.five_prime_matched_three_prime_missing;
        for (k, v) in &other.per_sample {
            *self.per_sample.entry(k.clone()).or_insert(0) += v;
        }
    }

    pub fn record_sample(&mut self, key: &SampleKey) {
        match key {
            SampleKey::Unassigned => {
                self.unassigned += 1;
            }
            other => {
                self.assigned += 1;
                *self.per_sample.entry(other.label()).or_insert(0) += 1;
            }
        }
    }
}

pub type RunStats = ChunkStats;

impl RunStats {
    pub fn print_summary(&self) {
        let t = self.total_reads.max(1) as f64;
        eprintln!();
        eprintln!("SeqMux summary");
        eprintln!("────────────────────────────────────────");
        eprintln!(
            "Total reads:        {:>12} ",
            format_count(self.total_reads)
        );
        eprintln!(
            "Assigned:           {:>12}  ({:5.2}%)",
            format_count(self.assigned),
            100.0 * self.assigned as f64 / t
        );
        eprintln!(
            "Unassigned:         {:>12}  ({:5.2}%)",
            format_count(self.unassigned),
            100.0 * self.unassigned as f64 / t
        );
        eprintln!(
            "Ambiguous:          {:>12}  ({:5.2}%)",
            format_count(self.ambiguous),
            100.0 * self.ambiguous as f64 / t
        );
        eprintln!(
            "Quality trimmed:    {:>12}  ({:5.2}%)",
            format_count(self.quality_trimmed),
            100.0 * self.quality_trimmed as f64 / t
        );
        eprintln!(
            "Adapter trimmed:    {:>12}  ({:5.2}%)",
            format_count(self.adapter_trimmed),
            100.0 * self.adapter_trimmed as f64 / t
        );
        eprintln!(
            "Length filtered:    {:>12}  ({:5.2}%)",
            format_count(self.too_short),
            100.0 * self.too_short as f64 / t
        );
        if self.five_prime_matched_three_prime_missing > 0 {
            eprintln!(
                "5' ok / 3' missing: {:>12}",
                format_count(self.five_prime_matched_three_prime_missing)
            );
        }
        if !self.per_sample.is_empty() {
            eprintln!();
            eprintln!("Per-sample counts:");
            for (name, count) in &self.per_sample {
                eprintln!("  {name:<40} {:>12}", format_count(*count));
            }
        }
        eprintln!();
    }

    pub fn write_tsv(&self, path: &Path) -> crate::error::Result<()> {
        let mut f = std::fs::File::create(path)?;
        let mut buf = String::new();
        writeln!(buf, "metric\tcount").ok();
        writeln!(buf, "total_reads\t{}", self.total_reads).ok();
        writeln!(buf, "assigned\t{}", self.assigned).ok();
        writeln!(buf, "unassigned\t{}", self.unassigned).ok();
        writeln!(buf, "ambiguous\t{}", self.ambiguous).ok();
        writeln!(buf, "quality_trimmed\t{}", self.quality_trimmed).ok();
        writeln!(buf, "adapter_trimmed\t{}", self.adapter_trimmed).ok();
        writeln!(buf, "too_short\t{}", self.too_short).ok();
        writeln!(
            buf,
            "five_prime_matched_three_prime_missing\t{}",
            self.five_prime_matched_three_prime_missing
        )
        .ok();
        for (name, count) in &self.per_sample {
            writeln!(buf, "sample:{name}\t{count}").ok();
        }
        f.write_all(buf.as_bytes())?;
        Ok(())
    }
}
