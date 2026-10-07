use seqmux::barcode::{assign_single, parse_barcodes_csv, MatchResult};

fn assigned(distance: usize) -> MatchResult {
    MatchResult::Match {
        sample_index: 0,
        distance,
        swapped: false,
    }
}

#[test]
fn single_barcode_margin_truth_table() {
    for (a, b, threshold, delta, expected) in [
        ("AAAA", "AACC", 3, 2, assigned(0)),
        ("AAAC", "AACC", 3, 2, MatchResult::Ambiguous { distance: 1 }),
        ("AAAC", "ACCC", 3, 2, assigned(1)),
        ("AAAC", "AACA", 3, 0, MatchResult::Ambiguous { distance: 1 }),
        ("AAAC", "AACC", 1, 0, assigned(1)), // delta=0 backward compat
        ("AAAC", "AACC", 1, 1, assigned(1)), // runner-up B has dist 2, delta=1 >= 1 -> assigned
        ("AAAC", "AACC", 1, 2, MatchResult::Ambiguous { distance: 1 }), // runner-up B (ineligible) triggers delta failure
        (
            "AAAC",
            "AACC",
            1,
            10,
            MatchResult::Ambiguous { distance: 1 },
        ),
        ("AAAC", "AACC", 0, 10, MatchResult::NoMatch),
        ("AAAC", "AACC", 3, 0, assigned(1)),
    ] {
        let csv = format!("SampleNumber,Barcode1\nA,{a}\nB,{b}\n");
        let mut cfg = parse_barcodes_csv(&csv, threshold, 0).unwrap();
        cfg.min_mismatch_delta = delta;
        assert_eq!(
            assign_single(b"AAAA", &cfg),
            expected,
            "{csv}, delta={delta}"
        );
    }
}

#[test]
fn ineligible_but_close_competitor() {
    // max1=1, max2=1, min_delta=2
    // Sample A: d1=1 d2=0 total=1 eligible
    // Sample B: d1=2 d2=0 total=2 ineligible
    // delta = 2 - 1 = 1 < 2 => Ambiguous
    let csv = "SampleNumber,Barcode1,Barcode2\nA,AAAC,CCCC\nB,AACC,CCCC\n";
    let mut cfg = parse_barcodes_csv(csv, 1, 1).unwrap();
    cfg.min_mismatch_delta = 2;
    // Read: AAAATTTTCCCC
    // For A: BC1=AAAC (dist 1 <= 1), BC2=CCCC (dist 0 <= 1) -> total 1, eligible
    // For B: BC1=AACC (dist 2 > 1, INELIGIBLE), BC2=CCCC (dist 0) -> total 2
    assert_eq!(
        assign_single(b"AAAATTTTCCCC", &cfg),
        MatchResult::Ambiguous { distance: 1 }
    );
}

#[test]
fn ineligible_and_sufficiently_distant_competitor() {
    // Sample A: d1=1 d2=0 total=1 eligible
    // Sample B: d1=4 d2=0 total=4 ineligible
    // min_delta=2 => delta = 4 - 1 = 3 >= 2 => Assigned A
    let csv = "SampleNumber,Barcode1,Barcode2\nA,AAAC,CCCC\nB,GGGG,CCCC\n";
    let mut cfg = parse_barcodes_csv(csv, 1, 1).unwrap();
    cfg.min_mismatch_delta = 2;
    assert_eq!(assign_single(b"AAAATTTTCCCC", &cfg), assigned(1));
}

#[test]
fn min_delta_zero_compatibility() {
    // The same close ineligible competitor must NOT change legacy behavior at delta=0
    let csv = "SampleNumber,Barcode1,Barcode2\nA,AAAC,CCCC\nB,AACC,CCCC\n";
    let mut cfg = parse_barcodes_csv(csv, 1, 1).unwrap();
    cfg.min_mismatch_delta = 0;
    assert_eq!(assign_single(b"AAAATTTTCCCC", &cfg), assigned(1));
}

#[test]
fn parser_defaults_to_backward_compatible_margin() {
    let cfg = parse_barcodes_csv("SampleNumber,Barcode1\nA,AAAC\nB,AACC\n", 2, 0).unwrap();
    assert_eq!(cfg.min_mismatch_delta, 0);
    assert_eq!(assign_single(b"AAAA", &cfg), assigned(1));
}

#[test]
fn dual_se_combines_eligible_barcode_distances() {
    let mut cfg = parse_barcodes_csv(
        "SampleNumber,Barcode1,Barcode2\nA,AAAC,CCCC\nB,AACA,CCCA\n",
        1,
        1,
    )
    .unwrap();
    for (delta, expected) in [
        (0, assigned(1)),
        (1, assigned(1)),
        (2, MatchResult::Ambiguous { distance: 1 }),
    ] {
        cfg.min_mismatch_delta = delta;
        assert_eq!(assign_single(b"AAAATTTTCCCC", &cfg), expected);
    }
}

#[test]
fn single_umi_bases_do_not_contribute_to_sample_margin() {
    let mut cfg = parse_barcodes_csv("SampleNumber,Barcode1\nA,ANAA\nB,CNCA\n", 2, 0).unwrap();
    for base in b"ACGTN" {
        let read = [b'A', *base, b'A', b'A'];
        cfg.min_mismatch_delta = 2;
        assert_eq!(assign_single(&read, &cfg), assigned(0));
        cfg.min_mismatch_delta = 3;
        assert_eq!(
            assign_single(&read, &cfg),
            MatchResult::Ambiguous { distance: 0 }
        );
    }
}

#[test]
fn dual_se_umi_bases_do_not_contribute_to_sample_margin() {
    let mut cfg = parse_barcodes_csv(
        "SampleNumber,Barcode1,Barcode2\nA,ANAA,CCNC\nB,CNCA,CANC\n",
        2,
        1,
    )
    .unwrap();
    cfg.min_mismatch_delta = 3;
    assert_eq!(assign_single(b"ATAACCTC", &cfg), assigned(0));
    cfg.min_mismatch_delta = 4;
    assert_eq!(
        assign_single(b"ATAACCTC", &cfg),
        MatchResult::Ambiguous { distance: 0 }
    );
}

#[test]
fn dual_se_equal_best_and_short_reads_remain_unassigned() {
    let mut cfg = parse_barcodes_csv(
        "SampleNumber,Barcode1,Barcode2\nA,AAAC,CCCC\nB,AAAA,CCCA\n",
        1,
        1,
    )
    .unwrap();
    for delta in [0, 2] {
        cfg.min_mismatch_delta = delta;
        assert_eq!(
            assign_single(b"AAAACCCC", &cfg),
            MatchResult::Ambiguous { distance: 1 }
        );
        assert_eq!(assign_single(b"AAA", &cfg), MatchResult::NoMatch);
    }
}
