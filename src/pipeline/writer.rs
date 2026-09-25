use super::worker::ProcessedChunk;
use crate::error::Result;
use crate::fastq::{check_output_path, FastqWriter};
use crate::output::{output_path_for_label, Mate, OutputKey, UNASSIGNED_SAMPLE_ID};
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
    sample_labels: Vec<String>,
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
        sample_labels: Vec<String>,
    ) -> Self {
        Self {
            out_dir: out_dir.to_path_buf(),
            prefix: prefix.to_string(),
            gzip,
            compression_level,
            force,
            sample_labels: sample_labels.clone(),
            writers: HashMap::new(),
            pending: BTreeMap::new(),
            next_expected: 0,
            stats: RunStats::new(sample_labels),
        }
    }

    pub fn push(&mut self, chunk: ProcessedChunk) -> Result<()> {
        if chunk.id < self.next_expected {
            return Err(crate::error::AppError::WorkerFailure(format!(
                "duplicate chunk id {}: already written (next expected: {})",
                chunk.id, self.next_expected
            )));
        }
        if self.pending.contains_key(&chunk.id) {
            return Err(crate::error::AppError::WorkerFailure(format!(
                "duplicate chunk id {}: already pending",
                chunk.id
            )));
        }
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
                let label = if key.sample_id == UNASSIGNED_SAMPLE_ID {
                    "unassigned"
                } else if let Some(l) = self.sample_labels.get(key.sample_id) {
                    l.as_str()
                } else {
                    "unknown"
                };
                let path =
                    output_path_for_label(&self.out_dir, &self.prefix, label, key.mate, self.gzip);
                check_output_path(&path, self.force)?;
                let w = FastqWriter::create(&path, self.gzip, self.compression_level)?;
                self.writers.insert(key, w);
            }
            let w = self.writers.get_mut(&key).unwrap();
            w.write_all(&data)?;
        }
        Ok(())
    }

    pub fn finish(mut self) -> Result<RunStats> {
        if !self.pending.is_empty() {
            let pending_ids: Vec<u64> = self.pending.keys().copied().collect();
            return Err(crate::error::AppError::WorkerFailure(format!(
                "ordered writer finished with {} pending chunks (ids: {:?}), but next expected was {}; missing chunk invariant violated",
                self.pending.len(),
                pending_ids,
                self.next_expected
            )));
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::ChunkStats;
    use tempfile::tempdir;

    fn make_chunk(id: u64) -> ProcessedChunk {
        ProcessedChunk {
            id,
            outputs: HashMap::new(),
            stats: ChunkStats::default(),
        }
    }

    #[test]
    fn ordered_writer_happy_path_out_of_order() {
        let dir = tempdir().unwrap();
        let mut writer = OrderedWriter::new(dir.path(), "test", false, 1, true, Vec::new());
        // Push 1, then 0, then 2
        assert!(writer.push(make_chunk(1)).is_ok());
        assert_eq!(writer.next_expected, 0);
        assert!(writer.push(make_chunk(0)).is_ok());
        assert_eq!(writer.next_expected, 2);
        assert!(writer.push(make_chunk(2)).is_ok());
        assert_eq!(writer.next_expected, 3);
        assert!(writer.finish().is_ok());
    }

    #[test]
    fn ordered_writer_rejects_duplicate_already_written() {
        let dir = tempdir().unwrap();
        let mut writer = OrderedWriter::new(dir.path(), "test", false, 1, true, Vec::new());
        assert!(writer.push(make_chunk(0)).is_ok());
        assert_eq!(writer.next_expected, 1);
        let err = writer.push(make_chunk(0)).unwrap_err();
        assert!(err.to_string().contains("already written"));
    }

    #[test]
    fn ordered_writer_rejects_duplicate_pending() {
        let dir = tempdir().unwrap();
        let mut writer = OrderedWriter::new(dir.path(), "test", false, 1, true, Vec::new());
        assert!(writer.push(make_chunk(2)).is_ok());
        let err = writer.push(make_chunk(2)).unwrap_err();
        assert!(err.to_string().contains("already pending"));
    }

    #[test]
    fn ordered_writer_rejects_missing_chunk_on_finish() {
        let dir = tempdir().unwrap();
        let mut writer = OrderedWriter::new(dir.path(), "test", false, 1, true, Vec::new());
        // Push 0 and 2; chunk 1 is missing
        assert!(writer.push(make_chunk(0)).is_ok());
        assert!(writer.push(make_chunk(2)).is_ok());
        let err = writer.finish().unwrap_err();
        assert!(err.to_string().contains("missing chunk invariant violated"));
    }
}
