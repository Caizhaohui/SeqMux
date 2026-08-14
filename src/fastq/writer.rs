use super::record::OwnedFastqRecord;
use crate::error::{AppError, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// Writer that can emit plain or gzip FASTQ.
pub enum FastqWriter {
    Plain(BufWriter<File>),
    Gzip(Box<BufWriter<GzEncoder<File>>>),
}

impl FastqWriter {
    pub fn create<P: AsRef<Path>>(path: P, gzip: bool, compression_level: u32) -> Result<Self> {
        let path = path.as_ref();
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)?;
            }
        }
        let file = File::create(path)?;
        if gzip {
            let level = Compression::new(compression_level.clamp(1, 9));
            let enc = GzEncoder::new(file, level);
            Ok(Self::Gzip(Box::new(BufWriter::with_capacity(
                256 * 1024,
                enc,
            ))))
        } else {
            Ok(Self::Plain(BufWriter::with_capacity(256 * 1024, file)))
        }
    }

    pub fn write_record(&mut self, rec: &OwnedFastqRecord) -> Result<()> {
        let bytes = rec.to_fastq_bytes();
        self.write_all(&bytes)
    }

    pub fn write_all(&mut self, data: &[u8]) -> Result<()> {
        match self {
            Self::Plain(w) => w.write_all(data)?,
            Self::Gzip(w) => w.write_all(data)?,
        }
        Ok(())
    }

    pub fn finish(self) -> Result<()> {
        match self {
            Self::Plain(mut w) => {
                w.flush()?;
            }
            Self::Gzip(w) => {
                let enc = (*w)
                    .into_inner()
                    .map_err(|e| AppError::Io(e.into_error()))?;
                enc.finish()?;
            }
        }
        Ok(())
    }
}

/// Check that path does not already exist unless force is set.
pub fn check_output_path(path: &Path, force: bool) -> Result<()> {
    if path.exists() && !force {
        return Err(AppError::OutputConflict {
            path: PathBuf::from(path),
        });
    }
    Ok(())
}
