use super::worker::ProcessedChunk;
use crate::error::Result;
use crate::fastq::{check_output_path, FastqWriter};
use crate::output::{output_path, Mate, OutputKey};
use crate::stats::RunStats;
use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

/// Ordered writer that preserves input order within each sample.
pub struct OrderedWriter {
    out_dir: PathBuf,
    prefix: String,
    gzip: bool,
    compression_level: u32,
    force: bool,
    writers: HashMap<OutputKey, FastqWriter>,
    pending: BTreeMap<u64, ProcessedChunk>,
    next_expected: u64,
    pub stats: RunStats,
}

impl OrderedWriter {
    pub fn new(
        out_dir: &Path,
        prefix: &str,
        gzip: bool,
        compression_level: u32,
        force: bool,
    ) -> Self {
        Self {
            out_dir: out_dir.to_path_buf(),
            prefix: prefix.to_string(),
            gzip,
            compression_level,
            force,
            writers: HashMap::new(),
            pending: BTreeMap::new(),
            next_expected: 0,
            stats: RunStats::default(),
        }
    }

    pub fn push(&mut self, chunk: ProcessedChunk) -> Result<()> {
        self.pending.insert(chunk.id, chunk);
        while let Some(chunk) = self.pending.remove(&self.next_expected) {
            self.write_chunk(chunk)?;
            self.next_expected += 1;
        }
        Ok(())
    }

    fn write_chunk(&mut self, chunk: ProcessedChunk) -> Result<()> {
        self.stats.merge(&chunk.stats);
        for (key, data) in chunk.outputs {
            if data.is_empty() {
                continue;
            }
            if !self.writers.contains_key(&key) {
                let path = output_path(
                    &self.out_dir,
                    &self.prefix,
                    &key.sample,
                    key.mate,
                    self.gzip,
                );
                check_output_path(&path, self.force)?;
                let w = FastqWriter::create(&path, self.gzip, self.compression_level)?;
                self.writers.insert(key.clone(), w);
            }
            let w = self.writers.get_mut(&key).unwrap();
            w.write_all(&data)?;
        }
        Ok(())
    }

    pub fn finish(mut self) -> Result<RunStats> {
        // Drain any remaining (should be empty if all chunks arrived)
        let remaining: Vec<u64> = self.pending.keys().copied().collect();
        for id in remaining {
            if let Some(chunk) = self.pending.remove(&id) {
                self.write_chunk(chunk)?;
            }
        }
        for (_, w) in self.writers.drain() {
            w.finish()?;
        }
        Ok(self.stats)
    }
}

/// Collect all sample keys that will be written (for pre-declaration if needed).
#[allow(dead_code)]
pub fn mate_for_paired(is_r1: bool) -> Mate {
    if is_r1 {
        Mate::R1
    } else {
        Mate::R2
    }
}
