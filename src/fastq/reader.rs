use super::record::{OwnedFastqRecord, ReadOrPair, ReadPair};
use crate::error::{AppError, Result};
use crate::util::normalize_read_id;
use crossbeam_channel::{bounded, Receiver, Sender};
use needletail::parse_fastx_file;
use needletail::FastxReader;
use std::collections::VecDeque;
use std::path::{Path, PathBuf};
use std::thread;

/// Trait for sequential FASTQ sources.
pub trait FastqSource {
    fn next_record(&mut self) -> Result<Option<OwnedFastqRecord>>;
}

/// Single-end FASTQ reader (plain or gzip).
pub struct FastqReader {
    reader: Box<dyn FastxReader>,
    index: u64,
}

impl FastqReader {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let reader = parse_fastx_file(path.as_ref()).map_err(|e| {
            AppError::FastqFormat(format!("failed to open {}: {e}", path.as_ref().display()))
        })?;
        Ok(Self { reader, index: 0 })
    }
}

impl FastqSource for FastqReader {
    fn next_record(&mut self) -> Result<Option<OwnedFastqRecord>> {
        match self.reader.next() {
            Some(Ok(rec)) => {
                self.index += 1;
                let name = rec.id().to_vec();
                let sequence: Vec<u8> = rec.seq().iter().map(|b| b.to_ascii_uppercase()).collect();
                let qualities = rec.qual().map(|q| q.to_vec()).ok_or_else(|| {
                    AppError::FastqFormat(format!(
                        "record {} missing quality scores: {}",
                        self.index,
                        String::from_utf8_lossy(&name)
                    ))
                })?;
                OwnedFastqRecord::new(name, sequence, qualities).map(Some)
            }
            Some(Err(e)) => Err(AppError::FastqFormat(format!(
                "parse error at record {}: {e}",
                self.index + 1
            ))),
            None => Ok(None),
        }
    }
}

/// Paired-end synchronized reader.
pub struct PairedFastqReader {
    r1: FastqReader,
    r2: FastqReader,
    index: u64,
}

impl PairedFastqReader {
    pub fn from_paths<P: AsRef<Path>, Q: AsRef<Path>>(r1: P, r2: Q) -> Result<Self> {
        Ok(Self {
            r1: FastqReader::from_path(r1)?,
            r2: FastqReader::from_path(r2)?,
            index: 0,
        })
    }

    pub fn next_pair(&mut self) -> Result<Option<ReadPair>> {
        let rec1 = self.r1.next_record()?;
        let rec2 = self.r2.next_record()?;
        match (rec1, rec2) {
            (None, None) => Ok(None),
            (Some(_), None) => Err(AppError::FastqFormat(format!(
                "R2 ended before R1 at record {}",
                self.index + 1
            ))),
            (None, Some(_)) => Err(AppError::FastqFormat(format!(
                "R1 ended before R2 at record {}",
                self.index + 1
            ))),
            (Some(r1), Some(r2)) => {
                self.index += 1;
                let id1 = normalize_read_id(&r1.name);
                let id2 = normalize_read_id(&r2.name);
                if id1 != id2 {
                    return Err(AppError::InvalidPair {
                        index: self.index,
                        r1: String::from_utf8_lossy(&r1.name).into_owned(),
                        r2: String::from_utf8_lossy(&r2.name).into_owned(),
                    });
                }
                Ok(Some(ReadPair { r1, r2 }))
            }
        }
    }
}

/// A chunk of FASTQ records read from a single input mate.
pub struct FastqRecordChunk {
    pub chunk_id: u64,
    pub records: Vec<OwnedFastqRecord>,
}

fn spawn_mate_reader_thread(
    path: PathBuf,
    chunk_reads: usize,
    skip_reads: u64,
    max_reads: u64,
    tx: Sender<Result<FastqRecordChunk>>,
) -> thread::JoinHandle<Result<u64>> {
    thread::spawn(move || {
        let run = || -> Result<u64> {
            let mut reader = FastqReader::from_path(&path)?;
            for _ in 0..skip_reads {
                if reader.next_record()?.is_none() {
                    return Ok(0);
                }
            }
            let mut chunk_id = 0u64;
            let mut n_read = 0u64;
            loop {
                let mut batch = Vec::with_capacity(chunk_reads);
                while batch.len() < chunk_reads {
                    if max_reads > 0 && n_read >= max_reads {
                        break;
                    }
                    match reader.next_record()? {
                        Some(rec) => {
                            n_read += 1;
                            batch.push(rec);
                        }
                        None => break,
                    }
                }
                if batch.is_empty() {
                    break;
                }
                if tx
                    .send(Ok(FastqRecordChunk {
                        chunk_id,
                        records: batch,
                    }))
                    .is_err()
                {
                    break;
                }
                chunk_id += 1;
                if max_reads > 0 && n_read >= max_reads {
                    break;
                }
            }
            Ok(chunk_id)
        };

        match run() {
            Ok(chunks) => Ok(chunks),
            Err(e) => {
                let _ = tx.send(Err(AppError::Other(e.to_string())));
                Err(e)
            }
        }
    })
}

/// Concurrent paired-end reader decompressing R1 and R2 concurrently on two separate threads.
pub struct ConcurrentPairedReader {
    r1_rx: Receiver<Result<FastqRecordChunk>>,
    r2_rx: Receiver<Result<FastqRecordChunk>>,
    r1_handle: Option<thread::JoinHandle<Result<u64>>>,
    r2_handle: Option<thread::JoinHandle<Result<u64>>>,
    next_expected_chunk_id: u64,
    total_pairs: u64,
    buffered_pairs: VecDeque<ReadPair>,
}

impl ConcurrentPairedReader {
    pub fn from_paths<P: AsRef<Path>, Q: AsRef<Path>>(
        r1: P,
        r2: Q,
        chunk_reads: usize,
        skip_reads: u64,
        max_reads: u64,
    ) -> Result<Self> {
        let p1 = r1.as_ref().to_path_buf();
        let p2 = r2.as_ref().to_path_buf();

        if !p1.exists() {
            return Err(AppError::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("input file not found: {}", p1.display()),
            )));
        }
        if !p2.exists() {
            return Err(AppError::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("input file not found: {}", p2.display()),
            )));
        }

        let (r1_tx, r1_rx) = bounded(8);
        let (r2_tx, r2_rx) = bounded(8);

        let r1_handle = spawn_mate_reader_thread(p1, chunk_reads, skip_reads, max_reads, r1_tx);
        let r2_handle = spawn_mate_reader_thread(p2, chunk_reads, skip_reads, max_reads, r2_tx);

        Ok(Self {
            r1_rx,
            r2_rx,
            r1_handle: Some(r1_handle),
            r2_handle: Some(r2_handle),
            next_expected_chunk_id: 0,
            total_pairs: 0,
            buffered_pairs: VecDeque::new(),
        })
    }

    /// Read the next matched chunk of paired reads.
    pub fn next_paired_chunk(&mut self) -> Result<Option<(u64, Vec<ReadOrPair>)>> {
        let c1 = self.r1_rx.recv();
        let c2 = self.r2_rx.recv();

        match (c1, c2) {
            (Ok(Err(e)), _) => Err(e),
            (_, Ok(Err(e))) => Err(e),
            (Ok(Ok(chunk1)), Ok(Ok(chunk2))) => {
                if chunk1.chunk_id != self.next_expected_chunk_id
                    || chunk2.chunk_id != self.next_expected_chunk_id
                {
                    return Err(AppError::WorkerFailure(format!(
                        "chunk ID mismatch: expected {}, got R1={}, R2={}",
                        self.next_expected_chunk_id, chunk1.chunk_id, chunk2.chunk_id
                    )));
                }
                if chunk1.records.len() != chunk2.records.len() {
                    return Err(AppError::FastqFormat(format!(
                        "R1 and R2 chunk record count mismatch at chunk {}: R1 has {}, R2 has {}",
                        chunk1.chunk_id,
                        chunk1.records.len(),
                        chunk2.records.len()
                    )));
                }
                let mut pairs = Vec::with_capacity(chunk1.records.len());
                for (r1, r2) in chunk1.records.into_iter().zip(chunk2.records) {
                    self.total_pairs += 1;
                    let id1 = normalize_read_id(&r1.name);
                    let id2 = normalize_read_id(&r2.name);
                    if id1 != id2 {
                        return Err(AppError::InvalidPair {
                            index: self.total_pairs,
                            r1: String::from_utf8_lossy(&r1.name).into_owned(),
                            r2: String::from_utf8_lossy(&r2.name).into_owned(),
                        });
                    }
                    pairs.push(ReadOrPair::Pair(ReadPair { r1, r2 }));
                }
                let id = chunk1.chunk_id;
                self.next_expected_chunk_id += 1;
                Ok(Some((id, pairs)))
            }
            (Err(_), Err(_)) => {
                self.finish()?;
                Ok(None)
            }
            (Ok(Ok(_)), Err(_)) => Err(AppError::FastqFormat(format!(
                "R2 ended before R1 at record {}",
                self.total_pairs + 1
            ))),
            (Err(_), Ok(Ok(_))) => Err(AppError::FastqFormat(format!(
                "R1 ended before R2 at record {}",
                self.total_pairs + 1
            ))),
        }
    }

    pub fn next_pair(&mut self) -> Result<Option<ReadPair>> {
        loop {
            if let Some(pair) = self.buffered_pairs.pop_front() {
                return Ok(Some(pair));
            }
            match self.next_paired_chunk()? {
                Some((_, items)) => {
                    for item in items {
                        if let ReadOrPair::Pair(p) = item {
                            self.buffered_pairs.push_back(p);
                        }
                    }
                    if !self.buffered_pairs.is_empty() {
                        return Ok(self.buffered_pairs.pop_front());
                    }
                }
                None => return Ok(None),
            }
        }
    }

    pub fn finish(&mut self) -> Result<()> {
        if let Some(h) = self.r1_handle.take() {
            h.join()
                .map_err(|_| AppError::WorkerFailure("R1 reader thread panicked".into()))??;
        }
        if let Some(h) = self.r2_handle.take() {
            h.join()
                .map_err(|_| AppError::WorkerFailure("R2 reader thread panicked".into()))??;
        }
        Ok(())
    }
}

/// Unified input iterator producing `ReadOrPair`.
pub enum InputReader {
    Single(FastqReader),
    Paired(PairedFastqReader),
    ConcurrentPaired(ConcurrentPairedReader),
}

impl InputReader {
    pub fn single<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(Self::Single(FastqReader::from_path(path)?))
    }

    pub fn paired<P: AsRef<Path>, Q: AsRef<Path>>(r1: P, r2: Q) -> Result<Self> {
        Ok(Self::Paired(PairedFastqReader::from_paths(r1, r2)?))
    }

    pub fn concurrent_paired<P: AsRef<Path>, Q: AsRef<Path>>(
        r1: P,
        r2: Q,
        chunk_reads: usize,
        skip_reads: u64,
        max_reads: u64,
    ) -> Result<Self> {
        Ok(Self::ConcurrentPaired(ConcurrentPairedReader::from_paths(
            r1,
            r2,
            chunk_reads,
            skip_reads,
            max_reads,
        )?))
    }

    pub fn next_item(&mut self) -> Result<Option<ReadOrPair>> {
        match self {
            Self::Single(r) => Ok(r.next_record()?.map(ReadOrPair::Single)),
            Self::Paired(r) => Ok(r.next_pair()?.map(ReadOrPair::Pair)),
            Self::ConcurrentPaired(r) => Ok(r.next_pair()?.map(ReadOrPair::Pair)),
        }
    }
}

/// Validate FASTQ file(s) and return summary statistics.
#[derive(Debug, Default)]
pub struct ValidateStats {
    pub records: u64,
    pub bases: u64,
    pub min_length: Option<usize>,
    pub max_length: Option<usize>,
    pub paired: bool,
}

pub fn validate_fastq<P: AsRef<Path>>(input: P, input2: Option<P>) -> Result<ValidateStats> {
    let mut stats = ValidateStats::default();
    if let Some(i2) = input2 {
        stats.paired = true;
        let mut reader = PairedFastqReader::from_paths(input, i2)?;
        while let Some(pair) = reader.next_pair()? {
            stats.records += 1;
            stats.bases += pair.r1.len() as u64 + pair.r2.len() as u64;
            for len in [pair.r1.len(), pair.r2.len()] {
                stats.min_length = Some(stats.min_length.map_or(len, |m| m.min(len)));
                stats.max_length = Some(stats.max_length.map_or(len, |m| m.max(len)));
            }
        }
    } else {
        let mut reader = FastqReader::from_path(input)?;
        while let Some(rec) = reader.next_record()? {
            stats.records += 1;
            stats.bases += rec.len() as u64;
            let len = rec.len();
            stats.min_length = Some(stats.min_length.map_or(len, |m| m.min(len)));
            stats.max_length = Some(stats.max_length.map_or(len, |m| m.max(len)));
        }
    }
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_fq(records: &[(&str, &str, &str)]) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        for (id, seq, qual) in records {
            writeln!(file, "@{id}\n{seq}\n+\n{qual}").unwrap();
        }
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_concurrent_paired_reader_happy_path() {
        let r1 = write_fq(&[
            ("read1/1", "ACGTACGT", "IIIIIIII"),
            ("read2/1", "TTTTAAAA", "JJJJJJJJ"),
            ("read3 1:N:0:1", "CCCCGGGG", "KKKKKKKK"),
        ]);
        let r2 = write_fq(&[
            ("read1/2", "TGCATGCA", "IIIIIIII"),
            ("read2/2", "AAAATTTT", "JJJJJJJJ"),
            ("read3 2:N:0:1", "GGGGCCCC", "KKKKKKKK"),
        ]);

        let mut reader = ConcurrentPairedReader::from_paths(r1.path(), r2.path(), 2, 0, 0).unwrap();
        let chunk1 = reader.next_paired_chunk().unwrap().expect("chunk 0");
        assert_eq!(chunk1.0, 0);
        assert_eq!(chunk1.1.len(), 2);

        let chunk2 = reader.next_paired_chunk().unwrap().expect("chunk 1");
        assert_eq!(chunk2.0, 1);
        assert_eq!(chunk2.1.len(), 1);

        assert!(reader.next_paired_chunk().unwrap().is_none());
        assert_eq!(reader.total_pairs, 3);
    }

    #[test]
    fn test_concurrent_paired_reader_id_mismatch() {
        let r1 = write_fq(&[
            ("read1/1", "ACGTACGT", "IIIIIIII"),
            ("read2/1", "TTTTAAAA", "JJJJJJJJ"),
        ]);
        let r2 = write_fq(&[
            ("read1/2", "TGCATGCA", "IIIIIIII"),
            ("read_mismatched/2", "AAAATTTT", "JJJJJJJJ"),
        ]);

        let mut reader =
            ConcurrentPairedReader::from_paths(r1.path(), r2.path(), 10, 0, 0).unwrap();
        let err = reader.next_paired_chunk().unwrap_err();
        match err {
            AppError::InvalidPair { index, .. } => assert_eq!(index, 2),
            other => panic!("expected InvalidPair, got {other:?}"),
        }
    }

    #[test]
    fn test_concurrent_paired_reader_r1_ended_early() {
        let r1 = write_fq(&[("read1/1", "ACGTACGT", "IIIIIIII")]);
        let r2 = write_fq(&[
            ("read1/2", "TGCATGCA", "IIIIIIII"),
            ("read2/2", "AAAATTTT", "JJJJJJJJ"),
        ]);

        let mut reader = ConcurrentPairedReader::from_paths(r1.path(), r2.path(), 1, 0, 0).unwrap();
        let _ = reader.next_paired_chunk().unwrap().unwrap();
        let err = reader.next_paired_chunk().unwrap_err();
        assert!(err.to_string().contains("R1 ended before R2"));
    }

    #[test]
    fn test_concurrent_paired_reader_r2_ended_early() {
        let r1 = write_fq(&[
            ("read1/1", "ACGTACGT", "IIIIIIII"),
            ("read2/1", "TTTTAAAA", "JJJJJJJJ"),
        ]);
        let r2 = write_fq(&[("read1/2", "TGCATGCA", "IIIIIIII")]);

        let mut reader = ConcurrentPairedReader::from_paths(r1.path(), r2.path(), 1, 0, 0).unwrap();
        let _ = reader.next_paired_chunk().unwrap().unwrap();
        let err = reader.next_paired_chunk().unwrap_err();
        assert!(err.to_string().contains("R2 ended before R1"));
    }

    #[test]
    fn test_concurrent_paired_reader_skip_and_max_reads() {
        let r1 = write_fq(&[
            ("read1/1", "ACGTACGT", "IIIIIIII"),
            ("read2/1", "TTTTAAAA", "JJJJJJJJ"),
            ("read3/1", "CCCCGGGG", "KKKKKKKK"),
            ("read4/1", "GGGGAAAA", "LLLLLLLL"),
        ]);
        let r2 = write_fq(&[
            ("read1/2", "TGCATGCA", "IIIIIIII"),
            ("read2/2", "AAAATTTT", "JJJJJJJJ"),
            ("read3/2", "GGGGCCCC", "KKKKKKKK"),
            ("read4/2", "CCCCGGGG", "LLLLLLLL"),
        ]);

        // skip 1, max 2 -> read2 and read3
        let mut reader =
            ConcurrentPairedReader::from_paths(r1.path(), r2.path(), 10, 1, 2).unwrap();
        let chunk = reader.next_paired_chunk().unwrap().unwrap();
        assert_eq!(chunk.1.len(), 2);
        if let ReadOrPair::Pair(p) = &chunk.1[0] {
            assert_eq!(p.r1.name, b"read2/1");
        }
        if let ReadOrPair::Pair(p) = &chunk.1[1] {
            assert_eq!(p.r1.name, b"read3/1");
        }
        assert!(reader.next_paired_chunk().unwrap().is_none());
    }
}
