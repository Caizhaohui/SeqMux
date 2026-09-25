/// Owned FASTQ record kept as bytes (no UTF-8 conversion in hot paths).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnedFastqRecord {
    pub name: Vec<u8>,
    pub sequence: Vec<u8>,
    pub qualities: Vec<u8>,
}

impl OwnedFastqRecord {
    pub fn new(name: Vec<u8>, sequence: Vec<u8>, qualities: Vec<u8>) -> crate::error::Result<Self> {
        if sequence.len() != qualities.len() {
            return Err(crate::error::AppError::FastqFormat(format!(
                "sequence length {} != quality length {} for read {}",
                sequence.len(),
                qualities.len(),
                String::from_utf8_lossy(&name)
            )));
        }
        Ok(Self {
            name,
            sequence,
            qualities,
        })
    }

    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Slice sequence and quality to [start, end).
    pub fn trim_to(&mut self, start: usize, end: usize) {
        let end = end.min(self.sequence.len());
        let start = start.min(end);
        self.sequence = self.sequence[start..end].to_vec();
        self.qualities = self.qualities[start..end].to_vec();
    }

    /// Drop the first `n` bases.
    pub fn trim_front(&mut self, n: usize) {
        let n = n.min(self.sequence.len());
        self.sequence.drain(..n);
        self.qualities.drain(..n);
    }

    /// Drop the last `n` bases.
    pub fn trim_back(&mut self, n: usize) {
        let n = n.min(self.sequence.len());
        let keep = self.sequence.len() - n;
        self.sequence.truncate(keep);
        self.qualities.truncate(keep);
    }

    /// Append UMI to header using Ultraplex-compatible `rbc:` tag.
    pub fn append_umi(&mut self, umi: &[u8]) {
        if umi.is_empty() {
            // Still add empty rbc: if not present? Ultraplex adds rbc: even for empty
            // when no match. For matched reads with UMI we always append.
            return;
        }
        // Replace spaces with underscores (Ultraplex behavior when adding UMI)
        for b in self.name.iter_mut() {
            if *b == b' ' {
                *b = b'_';
            }
        }
        if !self.name.windows(4).any(|w| w == b"rbc:") {
            self.name.extend_from_slice(b"rbc:");
        }
        self.name.extend_from_slice(umi);
    }

    /// Ensure `rbc:` tag exists (even if empty UMI) — matches Ultraplex no_match path.
    pub fn ensure_rbc_tag(&mut self) {
        if !self.name.windows(4).any(|w| w == b"rbc:") {
            self.name.extend_from_slice(b"rbc:");
        }
    }

    /// Append serialized FASTQ record directly to the provided buffer without temporary allocation.
    #[inline]
    pub fn append_fastq_to(&self, buf: &mut Vec<u8>) {
        buf.reserve(1 + self.name.len() + 1 + self.sequence.len() + 3 + self.qualities.len() + 1);
        buf.push(b'@');
        buf.extend_from_slice(&self.name);
        buf.push(b'\n');
        buf.extend_from_slice(&self.sequence);
        buf.extend_from_slice(b"\n+\n");
        buf.extend_from_slice(&self.qualities);
        buf.push(b'\n');
    }

    /// Serialize to FASTQ text (with trailing newline after quality).
    pub fn to_fastq_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(
            1 + self.name.len() + 1 + self.sequence.len() + 3 + self.qualities.len() + 1,
        );
        self.append_fastq_to(&mut buf);
        buf
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadPair {
    pub r1: OwnedFastqRecord,
    pub r2: OwnedFastqRecord,
}

#[derive(Debug, Clone)]
pub enum ReadOrPair {
    Single(OwnedFastqRecord),
    Pair(ReadPair),
}
