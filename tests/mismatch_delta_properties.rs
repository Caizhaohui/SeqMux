use proptest::prelude::*;
use seqmux::barcode::{
    assign_paired, assign_single, parse_barcodes_csv, BarcodeConfig, MatchResult, OrientationMode,
};
use std::collections::BTreeSet;

#[derive(Debug)]
struct Case {
    csv: String,
    r1: Vec<u8>,
    r2: Vec<u8>,
    m1: usize,
    m2: usize,
    path: usize,
    orientation: OrientationMode,
}

impl Case {
    fn config(&self) -> BarcodeConfig {
        parse_barcodes_csv(&self.csv, self.m1, self.m2).unwrap()
    }

    fn assign(&self, cfg: &BarcodeConfig) -> MatchResult {
        match self.path {
            0 | 1 => assign_single(&self.r1, cfg),
            2 => assign_paired(&self.r1, &self.r2, cfg, self.orientation),
            _ => unreachable!(),
        }
    }
}

fn cases() -> impl Strategy<Value = Case> {
    (
        prop::collection::btree_set(("[ACGT]{4}N?", "[ACGT]{4}N?"), 1..7),
        "[ACGTNacgt]{0,12}",
        "[ACGTNacgt]{0,12}",
        0usize..=4,
        0usize..=4,
        0usize..3,
        0usize..3,
    )
        .prop_map(|(rows, r1, r2, m1, m2, path, orientation)| {
            let mut csv = if path == 0 {
                "SampleNumber,Barcode1\n".to_string()
            } else {
                "SampleNumber,Barcode1,Barcode2\n".to_string()
            };
            let mut singles = BTreeSet::new();
            for (i, (b1, b2)) in rows.iter().enumerate() {
                if path == 0 {
                    if singles.insert(b1) {
                        csv.push_str(&format!("S{i},{b1}\n"));
                    }
                } else {
                    csv.push_str(&format!("S{i},{b1},{b2}\n"));
                }
            }
            Case {
                csv,
                r1: r1.into_bytes(),
                r2: r2.into_bytes(),
                m1,
                m2,
                path,
                orientation: [
                    OrientationMode::Both,
                    OrientationMode::Canonical,
                    OrientationMode::Swapped,
                ][orientation],
            }
        })
}

fn named(result: MatchResult, cfg: &BarcodeConfig) -> (Option<String>, Option<usize>, bool) {
    match result {
        MatchResult::Match {
            sample_index,
            distance,
            swapped,
        } => (
            Some(cfg.samples[sample_index].name.clone()),
            Some(distance),
            swapped,
        ),
        MatchResult::Ambiguous { distance } => (None, Some(distance), false),
        MatchResult::NoMatch => (None, None, false),
    }
}

fn raw_distance(read: &[u8], barcode: &[u8], suffix: bool) -> Option<usize> {
    if read.len() < barcode.len() {
        return None;
    }
    let start = if suffix {
        read.len() - barcode.len()
    } else {
        0
    };
    Some(
        read[start..start + barcode.len()]
            .iter()
            .zip(barcode)
            .filter(|(observed, expected)| {
                **expected != b'N' && observed.to_ascii_uppercase() != **expected
            })
            .count(),
    )
}

fn reference(case: &Case, cfg: &BarcodeConfig) -> MatchResult {
    // 1. enumerate every sample
    // 2. enumerate every allowed orientation
    // 3. calculate all distances
    // 4. determine each sample's best eligible orientation for assignment
    // 5. determine each sample's global minimum distance for confidence
    let mut sample_best_eligible: Vec<Option<(usize, bool)>> = vec![None; cfg.samples.len()];
    let mut sample_global_min: Vec<Option<usize>> = vec![None; cfg.samples.len()];

    for (i, sample) in cfg.samples.iter().enumerate() {
        let mut best_eligible: Option<(usize, bool)> = None;
        let mut min_global: Option<usize> = None;

        for swapped in [false, true] {
            if (case.path != 2 && swapped)
                || (case.path == 2
                    && matches!(
                        (case.orientation, swapped),
                        (OrientationMode::Canonical, true) | (OrientationMode::Swapped, false)
                    ))
            {
                continue;
            }
            let (r1, r2) = if swapped {
                (&case.r2, &case.r1)
            } else {
                (&case.r1, &case.r2)
            };
            let Some(d1) = raw_distance(r1, &sample.barcode1.raw, false) else {
                continue;
            };
            let d2 = match &sample.barcode2 {
                None => Some(0),
                Some(bc) => raw_distance(
                    if case.path == 1 { r1 } else { r2 },
                    &bc.raw,
                    case.path == 1,
                ),
            };
            if let Some(d2) = d2 {
                let total = d1 + d2;
                // Track global minimum for this sample across allowed orientations
                min_global = Some(min_global.map_or(total, |m| m.min(total)));

                // Check assignment eligibility
                if d1 <= cfg.mismatches_1 && d2 <= cfg.mismatches_2 {
                    match best_eligible {
                        None => best_eligible = Some((total, swapped)),
                        Some((best_d, _)) => {
                            if total < best_d {
                                best_eligible = Some((total, swapped));
                            }
                            // On tie, keep canonical (swapped == false)
                        }
                    }
                }
            }
        }

        sample_best_eligible[i] = best_eligible;
        sample_global_min[i] = min_global;
    }

    // 6. choose best eligible distinct sample
    let mut best_eligible_dist = usize::MAX;
    let mut best_sample_idx = usize::MAX;
    let mut best_swapped = false;
    let mut is_ambiguous = false;

    for (i, eligible) in sample_best_eligible.iter().enumerate() {
        if let Some((d, swapped)) = *eligible {
            if d < best_eligible_dist {
                best_eligible_dist = d;
                best_sample_idx = i;
                best_swapped = swapped;
                is_ambiguous = false;
            } else if d == best_eligible_dist {
                is_ambiguous = true;
            }
        }
    }

    if best_eligible_dist == usize::MAX {
        return MatchResult::NoMatch;
    }
    if is_ambiguous {
        return MatchResult::Ambiguous {
            distance: best_eligible_dist,
        };
    }

    // 7. find closest other distinct sample globally
    // 8. apply min-mismatch-delta only when delta > 0
    if cfg.min_mismatch_delta > 0 {
        let mut runner_up_dist = usize::MAX;
        for (i, global) in sample_global_min.iter().enumerate() {
            if i == best_sample_idx {
                continue; // distinct sample
            }
            if let Some(d) = *global {
                runner_up_dist = runner_up_dist.min(d);
            }
        }

        if runner_up_dist != usize::MAX
            && (runner_up_dist <= best_eligible_dist
                || runner_up_dist - best_eligible_dist < cfg.min_mismatch_delta)
        {
            return MatchResult::Ambiguous {
                distance: best_eligible_dist,
            };
        }
    }

    MatchResult::Match {
        sample_index: best_sample_idx,
        distance: best_eligible_dist,
        swapped: best_swapped,
    }
}

proptest! {
    #[test]
    fn default_and_positive_margins_match_independent_sorted_reference(case in cases()) {
        let mut cfg = case.config();
        prop_assert_eq!(cfg.min_mismatch_delta, 0);
        for delta in 0..=10 {
            cfg.min_mismatch_delta = delta;
            prop_assert_eq!(case.assign(&cfg), reference(&case, &cfg));
        }
    }

    #[test]
    fn sample_table_row_order_does_not_change_logical_outcome(case in cases(), delta in 0usize..=10) {
        let mut cfg = case.config();
        cfg.min_mismatch_delta = delta;
        let mut lines = case.csv.lines();
        let header = lines.next().unwrap();
        let mut rows: Vec<_> = lines.collect();
        rows.reverse();
        let mut reordered = parse_barcodes_csv(&format!("{header}\n{}\n", rows.join("\n")), case.m1, case.m2).unwrap();
        reordered.min_mismatch_delta = delta;
        prop_assert_eq!(named(case.assign(&cfg), &cfg), named(case.assign(&reordered), &reordered));
    }

    #[test]
    fn increasing_delta_cannot_create_assignments(case in cases()) {
        let mut cfg = case.config();
        let mut previous = case.assign(&cfg);
        for delta in 1..=10 {
            cfg.min_mismatch_delta = delta;
            let current = case.assign(&cfg);
            match &current {
                MatchResult::Match { .. } => prop_assert_eq!(&current, &previous),
                MatchResult::Ambiguous { .. } => prop_assert_ne!(&previous, &MatchResult::NoMatch),
                MatchResult::NoMatch => prop_assert_eq!(&previous, &MatchResult::NoMatch),
            }
            previous = current;
        }
    }

    #[test]
    fn unique_exact_packed_assignment_is_stable(
        b1 in "[ACGT]{8}", b2 in "[ACGT]{8}", delta in any::<usize>(),
    ) {
        let mut cfg = parse_barcodes_csv(&format!("SampleNumber,Barcode1,Barcode2\nA,{b1},{b2}\n"), 0, 0).unwrap();
        cfg.min_mismatch_delta = delta;
        prop_assert!(cfg.fast_exact_8bp_pe.is_some());
        prop_assert_eq!(assign_paired(b1.as_bytes(), b2.as_bytes(), &cfg, OrientationMode::Both),
            MatchResult::Match { sample_index: 0, distance: 0, swapped: false });
    }
}

#[test]
fn near_neighbor_row_order_preserves_name_or_ambiguity() {
    for rows in ["A,AAAC\nB,AACC\n", "B,AACC\nA,AAAC\n"] {
        let mut cfg = parse_barcodes_csv(&format!("SampleNumber,Barcode1\n{rows}"), 2, 0).unwrap();
        assert_eq!(
            named(assign_single(b"AAAA", &cfg), &cfg),
            (Some("A".into()), Some(1), false)
        );
        cfg.min_mismatch_delta = 2;
        assert_eq!(
            assign_single(b"AAAA", &cfg),
            MatchResult::Ambiguous { distance: 1 }
        );
    }
}
