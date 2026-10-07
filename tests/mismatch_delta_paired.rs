use seqmux::barcode::matcher::match_dual_paired;
use seqmux::barcode::{assign_paired, parse_barcodes_csv, MatchResult, OrientationMode};

fn assigned(distance: usize, swapped: bool) -> MatchResult {
    MatchResult::Match {
        sample_index: 0,
        distance,
        swapped,
    }
}

#[test]
fn paired_margin_uses_distinct_samples_after_orientation_selection() {
    for (rows, r1, r2, expected) in [
        ("A,AA,AC\nB,CC,CC\n", "AA", "AG", assigned(1, false)),
        ("A,AA,AC\nB,CC,CC\n", "AG", "AA", assigned(1, true)),
        ("A,AA,AC\nB,CC,CC\n", "AA", "AA", assigned(1, false)),
        (
            "A,AA,AC\nB,AA,AT\n",
            "AA",
            "AG",
            MatchResult::Ambiguous { distance: 1 },
        ),
        (
            "A,AA,AC\nB,AA,CC\n",
            "AA",
            "AG",
            MatchResult::Ambiguous { distance: 1 },
        ),
        ("A,AA,AC\n", "AA", "AG", assigned(1, false)),
    ] {
        let csv = format!("SampleNumber,Barcode1,Barcode2\n{rows}");
        let mut cfg = parse_barcodes_csv(&csv, 2, 2).unwrap();
        cfg.min_mismatch_delta = 2;
        assert!(cfg.fast_exact_8bp_pe.is_none());
        assert_eq!(
            assign_paired(r1.as_bytes(), r2.as_bytes(), &cfg, OrientationMode::Both),
            expected,
            "{csv} {r1}/{r2}"
        );
    }
}

#[test]
fn orientation_collapse_same_sample_never_competes_with_itself() {
    // Only one sample in table.
    // Canonical dist = 1, Swapped dist = 2.
    // The swapped orientation of Sample A must NEVER compete against Sample A as runner-up!
    let csv = "SampleNumber,Barcode1,Barcode2\nA,AAAC,CCCA\n";
    let mut cfg = parse_barcodes_csv(csv, 2, 2).unwrap();
    cfg.min_mismatch_delta = 5; // even with large delta, cannot be ambiguous against itself!
                                // r1: AAAA, r2: CCCA -> Canonical: d1=1, d2=0, total=1.
                                // Swapped: d1=3, d2=2, total=5.
    assert_eq!(
        assign_paired(b"AAAA", b"CCCA", &cfg, OrientationMode::Both),
        assigned(1, false)
    );
}

#[test]
fn distinct_sample_close_competitor_in_swapped_orientation() {
    // Sample A: canonical d1=1, d2=0, total=1 (eligible)
    // Sample B: canonical d1=4, d2=4, total=8; BUT swapped d1=1, d2=1, total=2 (ineligible or eligible)
    // Under OrientationMode::Both, Sample B's swapped orientation distance (2) must be considered!
    // With min_delta = 2, delta = 2 - 1 = 1 < 2 => Ambiguous
    let csv = "SampleNumber,Barcode1,Barcode2\nA,AAAA,CCCC\nB,GGGG,TTTT\n";
    let mut cfg = parse_barcodes_csv(csv, 1, 1).unwrap();
    cfg.min_mismatch_delta = 2;

    // Suppose read pair: r1=AAAC, r2=CCCC
    // For A canonical: BC1=AAAA (d1=1), BC2=CCCC (d2=0) -> total 1 (eligible)
    // Suppose sample B is defined such that swapped orientation is close:
    let csv2 = "SampleNumber,Barcode1,Barcode2\nA,AAAA,CCCC\nB,CCCC,AAAT\n";
    let mut cfg2 = parse_barcodes_csv(csv2, 1, 1).unwrap();
    cfg2.min_mismatch_delta = 2;
    // For read pair: r1=AAAC, r2=CCCC
    // For A: canonical BC1=AAAA (d1=1), BC2=CCCC (d2=0) -> total 1, eligible
    // For B: swapped (r2 against B's BC1, r1 against B's BC2):
    // r2(CCCC) vs BC1(CCCC) -> d1=0
    // r1(AAAC) vs BC2(AAAT) -> d2=1
    // Total swapped for B = 0 + 1 = 1!
    // Delta = 1 - 1 = 0 < 2 => Ambiguous!
    assert_eq!(
        assign_paired(b"AAAC", b"CCCC", &cfg2, OrientationMode::Both),
        MatchResult::Ambiguous { distance: 1 }
    );
}

#[test]
fn paired_eligibility_respects_thresholds_and_backward_compatibility() {
    for b in ["AACC,CCCC", "AAAA,CCAA"] {
        let csv = format!("SampleNumber,Barcode1,Barcode2\nA,AAAC,CCCA\nB,{b}\n");
        let mut cfg = parse_barcodes_csv(&csv, 1, 1).unwrap();
        // At delta=0: backward compatibility must assign Sample A
        cfg.min_mismatch_delta = 0;
        assert_eq!(
            assign_paired(b"AAAA", b"CCCC", &cfg, OrientationMode::Canonical),
            assigned(2, false)
        );
        assert_eq!(
            assign_paired(b"CCCC", b"AAAA", &cfg, OrientationMode::Swapped),
            assigned(2, true)
        );
        assert_eq!(
            assign_paired(b"AAAA", b"CCCC", &cfg, OrientationMode::Swapped),
            MatchResult::NoMatch
        );

        // At delta=10: ineligible B is close competitor (total dist 2), delta = 2 - 2 = 0 < 10 => Ambiguous!
        cfg.min_mismatch_delta = 10;
        assert_eq!(
            assign_paired(b"AAAA", b"CCCC", &cfg, OrientationMode::Canonical),
            MatchResult::Ambiguous { distance: 2 }
        );
    }
}

#[test]
fn paired_cross_sample_orientation_tie_is_ambiguous_even_at_delta_zero() {
    let mut cfg = parse_barcodes_csv(
        "SampleNumber,Barcode1,Barcode2\nA,AAAC,CCCC\nB,CCCA,AAAA\n",
        1,
        1,
    )
    .unwrap();
    for delta in [0, 2] {
        cfg.min_mismatch_delta = delta;
        assert_eq!(
            assign_paired(b"AAAA", b"CCCC", &cfg, OrientationMode::Both),
            MatchResult::Ambiguous { distance: 1 }
        );
        assert_eq!(
            assign_paired(b"AAA", b"CCCC", &cfg, OrientationMode::Both),
            MatchResult::NoMatch
        );
    }
}

#[test]
fn paired_umi_positions_are_ignored_in_both_orientations() {
    let mut cfg = parse_barcodes_csv(
        "SampleNumber,Barcode1,Barcode2\nA,ANAA,CCNC\nB,CNCA,CANC\n",
        2,
        1,
    )
    .unwrap();
    cfg.min_mismatch_delta = 3;
    assert_eq!(
        assign_paired(b"ATAA", b"CCTC", &cfg, OrientationMode::Both),
        assigned(0, false)
    );
    assert_eq!(
        assign_paired(b"CCTC", b"ATAA", &cfg, OrientationMode::Both),
        assigned(0, true)
    );
    cfg.min_mismatch_delta = 4;
    assert_eq!(
        assign_paired(b"ATAA", b"CCTC", &cfg, OrientationMode::Both),
        MatchResult::Ambiguous { distance: 0 }
    );
}

#[test]
fn exact_fast_path_preserves_unique_hits_and_cross_sample_ties_at_any_delta() {
    for (rows, r1, r2, expected) in [
        (
            "A,AAAAAAAA,CCCCCCCC\nB,GGGGGGGG,TTTTTTTT\n",
            b"AAAAAAAA",
            b"CCCCCCCC",
            assigned(0, false),
        ),
        (
            "A,AAAAAAAA,CCCCCCCC\nB,GGGGGGGG,TTTTTTTT\n",
            b"CCCCCCCC",
            b"AAAAAAAA",
            assigned(0, true),
        ),
        (
            "A,AAAAAAAA,AAAAAAAA\n",
            b"AAAAAAAA",
            b"AAAAAAAA",
            assigned(0, false),
        ),
        (
            "A,AAAAAAAA,CCCCCCCC\nB,CCCCCCCC,AAAAAAAA\n",
            b"AAAAAAAA",
            b"CCCCCCCC",
            MatchResult::Ambiguous { distance: 0 },
        ),
        (
            "A,AAAAAAAA,CCCCCCCC\n",
            b"AAAAAAAN",
            b"CCCCCCCC",
            MatchResult::NoMatch,
        ),
    ] {
        let csv = format!("SampleNumber,Barcode1,Barcode2\n{rows}");
        let mut cfg = parse_barcodes_csv(&csv, 0, 0).unwrap();
        // At delta=0, exact fast path is used
        cfg.min_mismatch_delta = 0;
        let fast = cfg
            .fast_exact_8bp_pe
            .as_ref()
            .expect("packed exact path retained at delta=0");
        assert_eq!(assign_paired(r1, r2, &cfg, OrientationMode::Both), expected);
        for orientation in [
            OrientationMode::Both,
            OrientationMode::Canonical,
            OrientationMode::Swapped,
        ] {
            assert_eq!(
                fast.match_paired(r1, r2, orientation),
                match_dual_paired(r1, r2, &cfg, orientation)
            );
        }
    }
}
