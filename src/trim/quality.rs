/// BWA-style quality trimming (same algorithm as Ultraplex / cutadapt `quality_trim_index`).
///
/// Returns `(start, end)` half-open bounds into the quality array.
/// Phred offset defaults to 33.
pub fn quality_trim_bounds(
    qualities: &[u8],
    cutoff_5: u8,
    cutoff_3: u8,
    phred_offset: u8,
) -> (usize, usize) {
    let len = qualities.len();
    if len == 0 {
        return (0, 0);
    }

    let mut start: usize = 0;
    if cutoff_5 > 0 {
        let mut s: i32 = 0;
        let mut max_qual: i32 = 0;
        for (i, &b) in qualities.iter().enumerate() {
            s += cutoff_5 as i32 - (b as i32 - phred_offset as i32);
            if s < 0 {
                break;
            }
            if s > max_qual {
                max_qual = s;
                start = i + 1;
            }
        }
    }

    let mut stop: usize = len;
    if cutoff_3 > 0 {
        let mut s: i32 = 0;
        let mut max_qual: i32 = 0;
        for i in (0..len).rev() {
            let q = qualities[i] as i32 - phred_offset as i32;
            s += cutoff_3 as i32 - q;
            if s < 0 {
                break;
            }
            if s > max_qual {
                max_qual = s;
                stop = i;
            }
        }
    }

    if start >= stop {
        return (0, 0);
    }
    (start, stop)
}

/// NextSeq-style trim: treat G as low quality from the 3' end.
pub fn nextseq_trim_index(
    sequence: &[u8],
    qualities: &[u8],
    cutoff: u8,
    phred_offset: u8,
) -> usize {
    let mut s: i32 = 0;
    let mut max_qual: i32 = 0;
    let mut max_i = qualities.len();
    for i in (0..qualities.len()).rev() {
        let mut q = qualities[i] as i32 - phred_offset as i32;
        if sequence[i].eq_ignore_ascii_case(&b'G') {
            q = cutoff as i32 - 1;
        }
        s += cutoff as i32 - q;
        if s < 0 {
            break;
        }
        if s > max_qual {
            max_qual = s;
            max_i = i;
        }
    }
    max_i
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_trim_high_quality() {
        let quals = b"IIIIIIII"; // Phred 40
        let (s, e) = quality_trim_bounds(quals, 20, 20, 33);
        assert_eq!((s, e), (0, 8));
    }

    #[test]
    fn trim_all_low() {
        let quals = b"!!!!!!!!"; // Phred 0
        let (s, e) = quality_trim_bounds(quals, 1, 1, 33);
        assert_eq!((s, e), (0, 0));
    }

    #[test]
    fn empty() {
        assert_eq!(quality_trim_bounds(b"", 10, 10, 33), (0, 0));
    }

    #[test]
    fn zero_cutoff() {
        let quals = b"!!!!!!!!";
        let (s, e) = quality_trim_bounds(quals, 0, 0, 33);
        assert_eq!((s, e), (0, 8));
    }

    #[test]
    fn trim_3p_only() {
        // High then low: IIII!!!!
        let quals = b"IIII!!!!";
        let (s, e) = quality_trim_bounds(quals, 0, 20, 33);
        assert_eq!(s, 0);
        assert!(e < 8);
        assert!(e >= 4);
    }
}
