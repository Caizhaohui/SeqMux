use crate::output::UNASSIGNED_SAMPLE_ID;
use crate::util::format_count;
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
    /// Dual-barcode PE assigned in Barcode1@R1 / Barcode2@R2 orientation.
    pub orientation_canonical: u64,
    /// Dual-barcode PE assigned in Barcode2@R1 / Barcode1@R2 orientation.
    pub orientation_swapped: u64,
    pub per_sample: Vec<u64>,
}

impl ChunkStats {
    pub fn new(sample_count: usize) -> Self {
        Self {
            per_sample: vec![0; sample_count],
            ..Default::default()
        }
    }

    pub fn merge(&mut self, other: &ChunkStats) {
        self.total_reads += other.total_reads;
        self.quality_trimmed += other.quality_trimmed;
        self.adapter_trimmed += other.adapter_trimmed;
        self.assigned += other.assigned;
        self.unassigned += other.unassigned;
        self.ambiguous += other.ambiguous;
        self.too_short += other.too_short;
        self.five_prime_matched_three_prime_missing += other.five_prime_matched_three_prime_missing;
        self.orientation_canonical += other.orientation_canonical;
        self.orientation_swapped += other.orientation_swapped;
        if self.per_sample.len() < other.per_sample.len() {
            self.per_sample.resize(other.per_sample.len(), 0);
        }
        for (i, &v) in other.per_sample.iter().enumerate() {
            self.per_sample[i] += v;
        }
    }

    #[inline(always)]
    pub fn record_sample(&mut self, sample_id: usize) {
        if sample_id == UNASSIGNED_SAMPLE_ID {
            self.unassigned += 1;
        } else {
            self.assigned += 1;
            if sample_id >= self.per_sample.len() {
                self.per_sample.resize(sample_id + 1, 0);
            }
            self.per_sample[sample_id] += 1;
        }
    }

    #[inline(always)]
    pub fn record_sample_opt(&mut self, sample_id: Option<usize>) {
        match sample_id {
            None => self.unassigned += 1,
            Some(idx) => {
                self.assigned += 1;
                if idx >= self.per_sample.len() {
                    self.per_sample.resize(idx + 1, 0);
                }
                self.per_sample[idx] += 1;
            }
        }
    }
}

#[derive(Debug, Default, Clone)]
pub struct RunStats {
    pub chunk_stats: ChunkStats,
    pub sample_labels: Vec<String>,
}

impl std::ops::Deref for RunStats {
    type Target = ChunkStats;
    fn deref(&self) -> &Self::Target {
        &self.chunk_stats
    }
}

impl std::ops::DerefMut for RunStats {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.chunk_stats
    }
}

impl RunStats {
    pub fn new(sample_labels: Vec<String>) -> Self {
        let sample_count = sample_labels.len();
        Self {
            chunk_stats: ChunkStats::new(sample_count),
            sample_labels,
        }
    }

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
        if self.orientation_canonical + self.orientation_swapped > 0 {
            eprintln!(
                "Orientation R1=BC1: {:>12}  ({:5.2}%)",
                format_count(self.orientation_canonical),
                100.0 * self.orientation_canonical as f64 / t
            );
            eprintln!(
                "Orientation R1=BC2: {:>12}  ({:5.2}%)",
                format_count(self.orientation_swapped),
                100.0 * self.orientation_swapped as f64 / t
            );
        }
        let mut sample_counts: Vec<(&str, u64)> = self
            .per_sample
            .iter()
            .enumerate()
            .filter(|(_, &c)| c > 0)
            .map(|(i, &c)| {
                let name = self
                    .sample_labels
                    .get(i)
                    .map(|s| s.as_str())
                    .unwrap_or("unknown");
                (name, c)
            })
            .collect();
        if !sample_counts.is_empty() {
            sample_counts.sort_by(|a, b| a.0.cmp(b.0));
            eprintln!();
            eprintln!("Per-sample counts:");
            for (name, count) in sample_counts {
                eprintln!("  {name:<40} {:>12}", format_count(count));
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
        writeln!(buf, "orientation_canonical\t{}", self.orientation_canonical).ok();
        writeln!(buf, "orientation_swapped\t{}", self.orientation_swapped).ok();
        let match_rate = if self.total_reads == 0 {
            0.0
        } else {
            self.assigned as f64 / self.total_reads as f64
        };
        writeln!(buf, "match_rate\t{match_rate:.6}").ok();
        let mut sample_counts: Vec<(&str, u64)> = self
            .per_sample
            .iter()
            .enumerate()
            .filter(|(_, &c)| c > 0)
            .map(|(i, &c)| {
                let name = self
                    .sample_labels
                    .get(i)
                    .map(|s| s.as_str())
                    .unwrap_or("unknown");
                (name, c)
            })
            .collect();
        sample_counts.sort_by(|a, b| a.0.cmp(b.0));
        for (name, count) in sample_counts {
            writeln!(buf, "sample:{name}\t{count}").ok();
        }
        f.write_all(buf.as_bytes())?;
        Ok(())
    }
}
