use super::config::{
    pack_8bp_2bit, BarcodeConfig, CompiledBarcode, DemuxMode, ExactDual8Matcher, OrientationMode,
    SampleEntry, SampleKey,
};

/// Match result for sample assignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MatchResult {
    Match {
        sample_index: usize,
        distance: usize,
        /// True when dual-barcode PE matched Barcode2@R1 + Barcode1@R2.
        swapped: bool,
    },
    Ambiguous {
        distance: usize,
    },
    NoMatch,
}

impl ExactDual8Matcher {
    /// Zero-allocation, O(1) table lookup for 8-bp dual-barcode PE matching.
    #[inline]
    pub fn match_paired(&self, r1: &[u8], r2: &[u8], orientation: OrientationMode) -> MatchResult {
        let (c1, c2) = match (pack_8bp_2bit(r1), pack_8bp_2bit(r2)) {
            (Some(a), Some(b)) => (a, b),
            _ => return MatchResult::NoMatch,
        };

        let allow_canonical = matches!(
            orientation,
            OrientationMode::Both | OrientationMode::Canonical
        );
        let allow_swapped = matches!(
            orientation,
            OrientationMode::Both | OrientationMode::Swapped
        );

        let hit_can = if allow_canonical {
            let key = ((c1 as u32) << 16) | (c2 as u32);
            self.table.get(&key).copied()
        } else {
            None
        };

        let hit_swap = if allow_swapped {
            let key = ((c2 as u32) << 16) | (c1 as u32);
            self.table.get(&key).copied()
        } else {
            None
        };

        match (hit_can, hit_swap) {
            (None, None) => MatchResult::NoMatch,
            (Some(idx), None) => MatchResult::Match {
                sample_index: idx,
                distance: 0,
                swapped: false,
            },
            (None, Some(idx)) => MatchResult::Match {
                sample_index: idx,
                distance: 0,
                swapped: true,
            },
            (Some(i1), Some(i2)) => {
                if i1 == i2 {
                    MatchResult::Match {
                        sample_index: i1,
                        distance: 0,
                        swapped: false,
                    }
                } else {
                    MatchResult::Ambiguous { distance: 0 }
                }
            }
        }
    }
}

/// Hamming distance between observed bases and expected (equal length).
pub fn hamming(observed: &[u8], expected: &[u8]) -> usize {
    observed
        .iter()
        .zip(expected.iter())
        .filter(|(a, b)| **a != **b)
        .count()
}

/// Extract bases at given positions from a sequence.
#[allow(dead_code)]
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

/// Distance of barcode against the 5' end of a read.
/// Zero-allocation in-place comparison.
#[inline]
pub fn barcode_distance_5p(seq: &[u8], bc: &CompiledBarcode) -> Option<usize> {
    if seq.len() < bc.pattern_len {
        return None;
    }
    let mut dist = 0;
    for (&pos, &expected) in bc.informative_positions.iter().zip(&bc.expected_bases) {
        if seq[pos].to_ascii_uppercase() != expected {
            dist += 1;
        }
    }
    Some(dist)
}

/// Distance of barcode against the 3' end of a read (pattern aligned to suffix).
/// Zero-allocation in-place comparison.
#[inline]
pub fn barcode_distance_3p(seq: &[u8], bc: &CompiledBarcode) -> Option<usize> {
    if seq.len() < bc.pattern_len {
        return None;
    }
    let start = seq.len() - bc.pattern_len;
    let mut dist = 0;
    for (&pos, &expected) in bc.informative_positions.iter().zip(&bc.expected_bases) {
        if seq[start + pos].to_ascii_uppercase() != expected {
            dist += 1;
        }
    }
    Some(dist)
}

/// Extract UMI bases from the 5' barcode region.
pub fn extract_umi_5p(seq: &[u8], bc: &CompiledBarcode) -> Vec<u8> {
    if bc.umi_positions.is_empty() {
        return Vec::new();
    }
    bc.umi_positions
        .iter()
        .filter_map(|&p| seq.get(p).copied())
        .map(|b| b.to_ascii_uppercase())
        .collect()
}

/// Extract UMI bases from the 3' barcode region.
pub fn extract_umi_3p(seq: &[u8], bc: &CompiledBarcode) -> Vec<u8> {
    if seq.len() < bc.pattern_len || bc.umi_positions.is_empty() {
        return Vec::new();
    }
    let start = seq.len() - bc.pattern_len;
    bc.umi_positions
        .iter()
        .filter_map(|&p| seq.get(start + p).copied())
        .map(|b| b.to_ascii_uppercase())
        .collect()
}

#[allow(dead_code)]
fn pick_best(candidates: &[(usize, usize)]) -> MatchResult {
    if candidates.is_empty() {
        return MatchResult::NoMatch;
    }
    let best = candidates.iter().map(|(_, d)| *d).min().unwrap();
    let winners: Vec<usize> = candidates
        .iter()
        .filter(|(_, d)| *d == best)
        .map(|(i, _)| *i)
        .collect();
    if winners.len() == 1 {
        MatchResult::Match {
            sample_index: winners[0],
            distance: best,
            swapped: false,
        }
    } else {
        MatchResult::Ambiguous { distance: best }
    }
}

#[allow(dead_code)]
fn pick_best_oriented(candidates: &[(usize, usize, bool)]) -> MatchResult {
    if candidates.is_empty() {
        return MatchResult::NoMatch;
    }
    let best = candidates.iter().map(|(_, d, _)| *d).min().unwrap();
    let at_best: Vec<(usize, bool)> = candidates
        .iter()
        .filter(|(_, d, _)| *d == best)
        .map(|(i, _, swapped)| (*i, *swapped))
        .collect();
    let mut unique_samples = Vec::new();
    for (i, _) in &at_best {
        if !unique_samples.contains(i) {
            unique_samples.push(*i);
        }
    }
    if unique_samples.len() != 1 {
        return MatchResult::Ambiguous { distance: best };
    }
    let sample_index = unique_samples[0];
    let swapped = at_best.iter().all(|(_, swapped)| *swapped);
    MatchResult::Match {
        sample_index,
        distance: best,
        swapped,
    }
}

/// Match single-end or R1-only using Barcode1 at 5'.
pub fn match_single(seq: &[u8], cfg: &BarcodeConfig) -> MatchResult {
    let mut best_dist = usize::MAX;
    let mut best_sample = usize::MAX;
    let mut is_ambiguous = false;

    for (i, sample) in cfg.samples.iter().enumerate() {
        if let Some(d) = barcode_distance_5p(seq, &sample.barcode1) {
            if d <= cfg.mismatches_1 {
                if d < best_dist {
                    best_dist = d;
                    best_sample = i;
                    is_ambiguous = false;
                } else if d == best_dist {
                    is_ambiguous = true;
                }
            }
        }
    }
    if best_dist == usize::MAX {
        MatchResult::NoMatch
    } else if is_ambiguous {
        MatchResult::Ambiguous {
            distance: best_dist,
        }
    } else {
        MatchResult::Match {
            sample_index: best_sample,
            distance: best_dist,
            swapped: false,
        }
    }
}

/// Match dual-barcode single-end: Barcode1 at 5', Barcode2 at 3'.
pub fn match_dual_single_end(seq: &[u8], cfg: &BarcodeConfig) -> MatchResult {
    let mut best_dist = usize::MAX;
    let mut best_sample = usize::MAX;
    let mut is_ambiguous = false;

    for (i, sample) in cfg.samples.iter().enumerate() {
        let Some(bc2) = sample.barcode2.as_ref() else {
            continue;
        };
        let Some(d1) = barcode_distance_5p(seq, &sample.barcode1) else {
            continue;
        };
        let Some(d2) = barcode_distance_3p(seq, bc2) else {
            continue;
        };
        if d1 <= cfg.mismatches_1 && d2 <= cfg.mismatches_2 {
            let d = d1 + d2;
            if d < best_dist {
                best_dist = d;
                best_sample = i;
                is_ambiguous = false;
            } else if d == best_dist {
                is_ambiguous = true;
            }
        }
    }
    if best_dist == usize::MAX {
        MatchResult::NoMatch
    } else if is_ambiguous {
        MatchResult::Ambiguous {
            distance: best_dist,
        }
    } else {
        MatchResult::Match {
            sample_index: best_sample,
            distance: best_dist,
            swapped: false,
        }
    }
}

/// Match dual-barcode paired-end fallback using running state (no heap allocations).
pub fn match_dual_paired(
    r1: &[u8],
    r2: &[u8],
    cfg: &BarcodeConfig,
    orientation: OrientationMode,
) -> MatchResult {
    let allow_canonical = matches!(
        orientation,
        OrientationMode::Both | OrientationMode::Canonical
    );
    let allow_swapped = matches!(
        orientation,
        OrientationMode::Both | OrientationMode::Swapped
    );

    let mut best_dist = usize::MAX;
    let mut best_sample = usize::MAX;
    let mut best_swapped = false;
    let mut is_ambiguous = false;

    for (i, sample) in cfg.samples.iter().enumerate() {
        let Some(bc2) = sample.barcode2.as_ref() else {
            continue;
        };
        if allow_canonical {
            if let (Some(d1), Some(d2)) = (
                barcode_distance_5p(r1, &sample.barcode1),
                barcode_distance_5p(r2, bc2),
            ) {
                if d1 <= cfg.mismatches_1 && d2 <= cfg.mismatches_2 {
                    let d = d1 + d2;
                    if d < best_dist {
                        best_dist = d;
                        best_sample = i;
                        best_swapped = false;
                        is_ambiguous = false;
                    } else if d == best_dist && best_sample != i {
                        is_ambiguous = true;
                    }
                }
            }
        }
        if allow_swapped {
            if let (Some(d1), Some(d2)) = (
                barcode_distance_5p(r2, &sample.barcode1),
                barcode_distance_5p(r1, bc2),
            ) {
                if d1 <= cfg.mismatches_1 && d2 <= cfg.mismatches_2 {
                    let d = d1 + d2;
                    if d < best_dist {
                        best_dist = d;
                        best_sample = i;
                        best_swapped = true;
                        is_ambiguous = false;
                    } else if d == best_dist && best_sample != i {
                        is_ambiguous = true;
                    }
                }
            }
        }
    }

    if best_dist == usize::MAX {
        MatchResult::NoMatch
    } else if is_ambiguous {
        MatchResult::Ambiguous {
            distance: best_dist,
        }
    } else {
        MatchResult::Match {
            sample_index: best_sample,
            distance: best_dist,
            swapped: best_swapped,
        }
    }
}

/// High-level assign for SE read.
pub fn assign_single(seq: &[u8], cfg: &BarcodeConfig) -> MatchResult {
    match cfg.mode {
        DemuxMode::SingleBarcode => match_single(seq, cfg),
        DemuxMode::DualBarcode => match_dual_single_end(seq, cfg),
    }
}

/// High-level assign for PE pair (routes to fast_exact_8bp_pe when available).
pub fn assign_paired(
    r1: &[u8],
    r2: &[u8],
    cfg: &BarcodeConfig,
    orientation: OrientationMode,
) -> MatchResult {
    if let Some(ref fast) = cfg.fast_exact_8bp_pe {
        return fast.match_paired(r1, r2, orientation);
    }
    match cfg.mode {
        DemuxMode::SingleBarcode => match_single(r1, cfg),
        DemuxMode::DualBarcode => match_dual_paired(r1, r2, cfg, orientation),
    }
}

pub fn sample_key(sample: &SampleEntry) -> SampleKey {
    SampleKey::Named(sample.name.clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::barcode::config::parse_barcodes_csv;

    #[test]
    fn dual_pe_exact() {
        let csv = "SampleNumber,Barcode1,Barcode2\n\
                   A1,AAGTCCAA,GGAGTACT\n\
                   A2,GACCTGAA,GGACTTGG\n";
        let cfg = parse_barcodes_csv(csv, 0, 0).unwrap();
        let r1 = b"AAGTCCAAAAAAAAAAA";
        let r2 = b"GGAGTACTCCCCCCCC";
        let m = assign_paired(r1, r2, &cfg, OrientationMode::Both);
        assert_eq!(
            m,
            MatchResult::Match {
                sample_index: 0,
                distance: 0,
                swapped: false,
            }
        );
    }

    #[test]
    fn dual_pe_swapped_orientation() {
        let csv = "SampleNumber,Barcode1,Barcode2\n\
                   A1,AAGTCCAA,GGAGTACT\n\
                   A2,GACCTGAA,GGACTTGG\n";
        let cfg = parse_barcodes_csv(csv, 0, 0).unwrap();
        // mates swapped vs sample table
        let r1 = b"GGAGTACTCCCCCCCC";
        let r2 = b"AAGTCCAAAAAAAAAAA";
        let m = assign_paired(r1, r2, &cfg, OrientationMode::Both);
        assert_eq!(
            m,
            MatchResult::Match {
                sample_index: 0,
                distance: 0,
                swapped: true,
            }
        );
        let canonical_only = assign_paired(r1, r2, &cfg, OrientationMode::Canonical);
        assert_eq!(canonical_only, MatchResult::NoMatch);
    }

    #[test]
    fn dual_pe_mismatch() {
        let csv = "SampleNumber,Barcode1,Barcode2\ns1,AAAAAAAA,TTTTTTTT\n";
        let cfg = parse_barcodes_csv(csv, 1, 0).unwrap();
        // one mismatch in barcode1
        let r1 = b"AAAAAAATCCCCCCCC";
        let r2 = b"TTTTTTTTGGGGGGGG";
        let m = assign_paired(r1, r2, &cfg, OrientationMode::Both);
        assert!(matches!(m, MatchResult::Match { .. }));
    }

    #[test]
    fn single_barcode() {
        let csv = "SampleNumber,Barcode1\ns1,ATGATGAT\ns2,CCGTAACG\n";
        let cfg = parse_barcodes_csv(csv, 0, 0).unwrap();
        let m = assign_single(b"ATGATGATAAAAAAAA", &cfg);
        assert_eq!(
            m,
            MatchResult::Match {
                sample_index: 0,
                distance: 0,
                swapped: false,
            }
        );
    }

    #[test]
    fn no_match() {
        let csv = "SampleNumber,Barcode1,Barcode2\ns1,AAAAAAAA,TTTTTTTT\n";
        let cfg = parse_barcodes_csv(csv, 0, 0).unwrap();
        let m = assign_paired(
            b"GGGGGGGGCCCCCCCC",
            b"CCCCCCCCGGGGGGGG",
            &cfg,
            OrientationMode::Both,
        );
        assert_eq!(m, MatchResult::NoMatch);
    }
}
