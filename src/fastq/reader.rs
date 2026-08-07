use super::record::{OwnedFastqRecord, ReadOrPair, ReadPair};
use crate::error::{AppError, Result};
use crate::util::normalize_read_id;
use needletail::parse_fastx_file;
use needletail::FastxReader;
use std::path::Path;

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

/// Unified input iterator producing `ReadOrPair`.
pub enum InputReader {
    Single(FastqReader),
    Paired(PairedFastqReader),
}

impl InputReader {
    pub fn single<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(Self::Single(FastqReader::from_path(path)?))
    }

    pub fn paired<P: AsRef<Path>, Q: AsRef<Path>>(r1: P, r2: Q) -> Result<Self> {
        Ok(Self::Paired(PairedFastqReader::from_paths(r1, r2)?))
    }

    pub fn next_item(&mut self) -> Result<Option<ReadOrPair>> {
        match self {
            Self::Single(r) => Ok(r.next_record()?.map(ReadOrPair::Single)),
            Self::Paired(r) => Ok(r.next_pair()?.map(ReadOrPair::Pair)),
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
