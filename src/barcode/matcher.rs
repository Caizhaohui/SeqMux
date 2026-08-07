use super::config::{CompiledFivePrime, CompiledThreePrime};

/// Match result for a barcode search.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MatchResult {
    Match { index: usize, distance: usize },
    Ambiguous { distance: usize },
    NoMatch,
}

/// Hamming distance between observed bases and expected (equal length).
pub fn hamming(observed: &[u8], expected: &[u8]) -> usize {
    observed
        .iter()
        .zip(expected.iter())
        .filter(|(a, b)| {
            // N in the read counts as mismatch (Ultraplex: penalty for N in the read)
            **a != **b
        })
        .count()
}

/// Extract bases at given positions from a sequence.
pub fn extract_at(seq: &[u8], positions: &[usize]) -> Option<Vec<u8>> {
    let mut out = Vec::with_capacity(positions.len());
    for &p in positions {
        if p >= seq.len() {
            return None;
        }
        out.push(seq[p].to_ascii_uppercase());
    }
    Some(out)
}

/// Extract bases using offsets from the 3' end (0 = last base).
/// Returns bases in 5'→3' pattern order (same order as `expected_bases`).
pub fn extract_from_end(
    seq: &[u8],
    offsets_from_end: &[usize],
    pattern_order_indices: &[usize],
) -> Option<Vec<u8>> {
    // pattern_order_indices not needed if offsets are stored in pattern 5'→3' order.
    // Our compile_three stores informative_offsets_from_end in pattern 5'→3' order.
    let _ = pattern_order_indices;
    let mut out = Vec::with_capacity(offsets_from_end.len());
    for &off in offsets_from_end {
        if off >= seq.len() {
            return None;
        }
        let idx = seq.len() - 1 - off;
        out.push(seq[idx].to_ascii_uppercase());
    }
    Some(out)
}

/// Best match among 5' barcodes.
pub fn match_five_prime(
    seq: &[u8],
    barcodes: &[CompiledFivePrime],
    max_mismatches: usize,
) -> MatchResult {
    if barcodes.is_empty() {
        return MatchResult::NoMatch;
    }
    // All share same informative positions
    let positions = &barcodes[0].informative_positions;
    let observed = match extract_at(seq, positions) {
        Some(o) => o,
        None => return MatchResult::NoMatch,
    };

    let mut best_distance = usize::MAX;
    let mut best_idx: Option<usize> = None;
    let mut tie = false;

    for (i, bc) in barcodes.iter().enumerate() {
        let d = hamming(&observed, &bc.expected_bases);
        if d < best_distance {
            best_distance = d;
            best_idx = Some(i);
            tie = false;
        } else if d == best_distance {
            tie = true;
        }
    }

    if best_distance > max_mismatches {
        MatchResult::NoMatch
    } else if tie {
        MatchResult::Ambiguous {
            distance: best_distance,
        }
    } else {
        MatchResult::Match {
            index: best_idx.unwrap(),
            distance: best_distance,
        }
    }
}

/// Best match among linked 3' barcodes for one 5' group.
pub fn match_three_prime(
    seq: &[u8],
    barcodes: &[CompiledThreePrime],
    max_mismatches: usize,
) -> MatchResult {
    if barcodes.is_empty() {
        return MatchResult::NoMatch;
    }

    let mut best_distance = usize::MAX;
    let mut best_idx: Option<usize> = None;
    let mut tie = false;

    for (i, bc) in barcodes.iter().enumerate() {
        if seq.len() < bc.pattern_len {
            continue;
        }
        let observed = match extract_from_end(seq, &bc.informative_offsets_from_end, &[]) {
            Some(o) => o,
            None => continue,
        };
        let d = hamming(&observed, &bc.expected_bases);
        if d < best_distance {
            best_distance = d;
            best_idx = Some(i);
            tie = false;
        } else if d == best_distance {
            tie = true;
        }
    }

    match best_idx {
        Some(index) if best_distance <= max_mismatches && !tie => MatchResult::Match {
            index,
            distance: best_distance,
        },
        Some(_) if best_distance <= max_mismatches && tie => MatchResult::Ambiguous {
            distance: best_distance,
        },
        _ => MatchResult::NoMatch,
    }
}

/// Extract 5' UMI bases from sequence given compiled barcode.
pub fn extract_five_umi(seq: &[u8], bc: &CompiledFivePrime) -> Vec<u8> {
    bc.umi_positions
        .iter()
        .filter_map(|&p| seq.get(p).copied())
        .map(|b| b.to_ascii_uppercase())
        .collect()
}

/// Extract 3' UMI bases from sequence given compiled barcode.
pub fn extract_three_umi(seq: &[u8], bc: &CompiledThreePrime) -> Vec<u8> {
    // Return UMI in 5'→3' order along the pattern
    // umi_offsets_from_end are stored in pattern 5'→3' order
    let mut out = Vec::with_capacity(bc.umi_offsets_from_end.len());
    for &off in &bc.umi_offsets_from_end {
        if off < seq.len() {
            let idx = seq.len() - 1 - off;
            out.push(seq[idx].to_ascii_uppercase());
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::barcode::config::parse_barcodes_csv;

    #[test]
    fn exact_match() {
        let cfg = parse_barcodes_csv("NNNATGNN:a\nNNNCCGNN:b\n", 0, 0, false).unwrap();
        let seq = b"ACGATGTCAAAAAAAA";
        let r = match_five_prime(seq, &cfg.five_prime, 0);
        assert_eq!(
            r,
            MatchResult::Match {
                index: 0,
                distance: 0
            }
        );
        let umi = extract_five_umi(seq, &cfg.five_prime[0]);
        assert_eq!(umi, b"ACGTC");
    }

    #[test]
    fn mismatch_and_tie() {
        let cfg = parse_barcodes_csv("ATG:a\nATC:b\n", 1, 0, false).unwrap();
        // ATA is distance 1 from both ATG and ATC → ambiguous
        let r = match_five_prime(b"ATAAAA", &cfg.five_prime, 1);
        assert!(matches!(r, MatchResult::Ambiguous { .. }));
    }

    #[test]
    fn no_match_beyond_threshold() {
        let cfg = parse_barcodes_csv("ATG:a\n", 0, 0, false).unwrap();
        let r = match_five_prime(b"GGGAAAA", &cfg.five_prime, 0);
        assert_eq!(r, MatchResult::NoMatch);
    }

    #[test]
    fn three_prime_match() {
        let cfg = parse_barcodes_csv("ATG,NNTTCNN:s1\n", 0, 0, false).unwrap();
        let three = &cfg.five_prime[0].linked_three_prime;
        // pattern NNTTCNN, length 7; informative TTC at offsets from end
        // want ...NNTTCNN at end: e.g. AATTCGG
        let seq = b"XXXXXXAATTCGG";
        let r = match_three_prime(seq, three, 0);
        assert_eq!(
            r,
            MatchResult::Match {
                index: 0,
                distance: 0
            }
        );
    }
}
