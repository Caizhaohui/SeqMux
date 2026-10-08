use assert_cmd::Command;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use tempfile::tempdir;

fn write_fastq(path: &std::path::Path, records: &[(&str, &str, &str)]) {
    let mut f = fs::File::create(path).unwrap();
    for (name, seq, qual) in records {
        writeln!(f, "@{name}").unwrap();
        writeln!(f, "{seq}").unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "{qual}").unwrap();
    }
}

fn parse_summary_tsv(tsv_content: &str) -> (Vec<String>, HashMap<String, String>) {
    let mut lines = Vec::new();
    let mut map = HashMap::new();
    for line in tsv_content.lines() {
        lines.push(line.to_string());
        if let Some((k, v)) = line.split_once('\t') {
            map.insert(k.to_string(), v.to_string());
        }
    }
    (lines, map)
}

fn check_invariants(map: &HashMap<String, String>, is_dual_pe: bool) {
    let total_reads: u64 = map
        .get("total_reads")
        .expect("total_reads")
        .parse()
        .unwrap();
    let assigned: u64 = map.get("assigned").expect("assigned").parse().unwrap();
    let unassigned: u64 = map.get("unassigned").expect("unassigned").parse().unwrap();
    let ambiguous: u64 = map.get("ambiguous").expect("ambiguous").parse().unwrap();
    let too_short: u64 = map.get("too_short").expect("too_short").parse().unwrap();
    let matched_before_filter: u64 = map
        .get("matched_before_filter")
        .expect("matched_before_filter")
        .parse()
        .unwrap();
    let no_match_before_filter: u64 = map
        .get("no_match_before_filter")
        .expect("no_match_before_filter")
        .parse()
        .unwrap();
    let ambiguous_after_filter: u64 = map
        .get("ambiguous_after_filter")
        .expect("ambiguous_after_filter")
        .parse()
        .unwrap();
    let orientation_canonical: u64 = map
        .get("orientation_canonical")
        .expect("orientation_canonical")
        .parse()
        .unwrap();
    let orientation_swapped: u64 = map
        .get("orientation_swapped")
        .expect("orientation_swapped")
        .parse()
        .unwrap();

    let sum_samples: u64 = map
        .iter()
        .filter(|(k, _)| k.starts_with("sample:"))
        .map(|(_, v)| v.parse::<u64>().unwrap())
        .sum();

    // Invariant 1: total_reads == matched_before_filter + no_match_before_filter + ambiguous
    assert_eq!(
        total_reads,
        matched_before_filter + no_match_before_filter + ambiguous,
        "Invariant 1 failed: total_reads != matched_before_filter + no_match_before_filter + ambiguous"
    );

    // Invariant 2: total_reads == assigned + unassigned + too_short
    assert_eq!(
        total_reads,
        assigned + unassigned + too_short,
        "Invariant 2 failed: total_reads != assigned + unassigned + too_short"
    );

    // Invariant 3: sum(per_sample) == assigned
    assert_eq!(
        sum_samples, assigned,
        "Invariant 3 failed: sum(per_sample) != assigned"
    );

    // Invariant 4: ambiguous_after_filter <= unassigned
    assert!(
        ambiguous_after_filter <= unassigned,
        "Invariant 4 failed: ambiguous_after_filter ({ambiguous_after_filter}) > unassigned ({unassigned})"
    );

    // Invariant 5: when too_short == 0:
    if too_short == 0 {
        assert_eq!(
            matched_before_filter, assigned,
            "Invariant 5a failed: matched_before_filter != assigned when too_short == 0"
        );
        assert_eq!(
            no_match_before_filter + ambiguous,
            unassigned,
            "Invariant 5b failed: no_match_before_filter + ambiguous != unassigned when too_short == 0"
        );
        assert_eq!(
            ambiguous_after_filter, ambiguous,
            "Invariant 5c failed: ambiguous_after_filter != ambiguous when too_short == 0"
        );
    }

    // Invariant 6: Dual PE: orientation_canonical + orientation_swapped == matched_before_filter
    if is_dual_pe {
        assert_eq!(
            orientation_canonical + orientation_swapped,
            matched_before_filter,
            "Invariant 6 failed: canonical + swapped != matched_before_filter in dual PE"
        );
    }
}

#[test]
fn test_stats_invariants_single_end_matrix() {
    let dir = tempdir().unwrap();
    let fq = dir.path().join("se.fastq");
    let bc = dir.path().join("barcodes.csv");

    // Barcodes: 8bp
    // sample1: AAAAAAAA
    // sample2: TTTTTTTT
    // sample3: CCCCCCCC
    fs::write(
        &bc,
        "SampleNumber,Barcode1\n\
         sample1,AAAAAAAA\n\
         sample2,TTTTTTTT\n\
         sample3,CCCCCCCC\n",
    )
    .unwrap();

    // Construct 6 reads:
    // 1. Long unique match to sample1 (barcode AAAAAAAA, trimmed len 12 >= 10)
    // 2. Short unique match to sample1 (barcode AAAAAAAA, trimmed len 4 < 10)
    // 3. Long ambiguous (equal dist to AAAAAAAA and TTTTTTTT: e.g. AAAATTTT, trimmed len 12 >= 10)
    // 4. Short ambiguous (equal dist to AAAAAAAA and TTTTTTTT: e.g. AAAATTTT, trimmed len 4 < 10)
    // 5. Long no-match (GGGGGGGG, trimmed len 12 >= 10)
    // 6. Short no-match (GGGGGGGG, trimmed len 4 < 10)
    write_fastq(
        &fq,
        &[
            (
                "read1_long_match",
                "AAAAAAAACCCCCCCCCCCC",
                "IIIIIIIIIIIIIIIIIIII",
            ),
            ("read2_short_match", "AAAAAAAAACCC", "IIIIIIIIIIII"),
            (
                "read3_long_ambig",
                "AAAATTTTCCCCCCCCCCCC",
                "IIIIIIIIIIIIIIIIIIII",
            ),
            ("read4_short_ambig", "AAAATTTT", "IIIIIIII"),
            (
                "read5_long_nomatch",
                "GGGGGGGGCCCCCCCCCCCC",
                "IIIIIIIIIIIIIIIIIIII",
            ),
            ("read6_short_nomatch", "GGGGGGGG", "IIIIIIII"),
        ],
    );

    // Test combinations: min_length in [0, 10], threads in [1, 2], counts_only in [false, true], discard_unassigned in [false, true]
    for min_length in [0, 10] {
        for threads in [1, 2] {
            for counts_only in [false, true] {
                for discard_unassigned in [false, true] {
                    let out_dir = dir.path().join(format!(
                        "out_se_l{min_length}_t{threads}_c{counts_only}_d{discard_unassigned}"
                    ));

                    let mut cmd = Command::cargo_bin("seqmux").unwrap();
                    cmd.args([
                        "demux",
                        "-i",
                        fq.to_str().unwrap(),
                        "-b",
                        bc.to_str().unwrap(),
                        "-o",
                        out_dir.to_str().unwrap(),
                        "--min-length",
                        &min_length.to_string(),
                        "--threads",
                        &threads.to_string(),
                        "--mismatches-1",
                        "4",
                        "--no-adapter",
                        "--no-gzip",
                    ]);
                    if counts_only {
                        cmd.arg("--counts-only");
                    }
                    if discard_unassigned {
                        cmd.arg("--discard-unassigned");
                    }

                    cmd.assert().success();

                    let summary_file = out_dir.join("seqmux.summary.tsv");
                    let content = fs::read_to_string(&summary_file).unwrap();
                    let (lines, map) = parse_summary_tsv(&content);

                    check_invariants(&map, false);

                    // Check exact expected counts
                    let total: u64 = map["total_reads"].parse().unwrap();
                    assert_eq!(total, 6);

                    let matched_before: u64 = map["matched_before_filter"].parse().unwrap();
                    let no_match_before: u64 = map["no_match_before_filter"].parse().unwrap();
                    let ambig_before: u64 = map["ambiguous"].parse().unwrap();
                    assert_eq!(matched_before, 2);
                    assert_eq!(ambig_before, 2);
                    assert_eq!(no_match_before, 2);

                    let too_short: u64 = map["too_short"].parse().unwrap();
                    let assigned: u64 = map["assigned"].parse().unwrap();
                    let unassigned: u64 = map["unassigned"].parse().unwrap();
                    let ambig_after: u64 = map["ambiguous_after_filter"].parse().unwrap();

                    if min_length == 0 {
                        assert_eq!(too_short, 0);
                        assert_eq!(assigned, 2);
                        assert_eq!(unassigned, 4);
                        assert_eq!(ambig_after, 2);
                    } else {
                        // min_length == 10: 3 short reads filtered out (1 match, 1 ambig, 1 no-match)
                        assert_eq!(too_short, 3);
                        assert_eq!(assigned, 1);
                        assert_eq!(unassigned, 2);
                        assert_eq!(ambig_after, 1);
                    }

                    // Check TSV suffix positioning: last 4 lines must be the new metrics
                    let n = lines.len();
                    assert!(lines[n - 4].starts_with("matched_before_filter\t"));
                    assert!(lines[n - 3].starts_with("no_match_before_filter\t"));
                    assert!(lines[n - 2].starts_with("ambiguous_after_filter\t"));
                    assert!(lines[n - 1].starts_with("assignment_rate_before_filter\t"));
                }
            }
        }
    }
}

#[test]
fn test_stats_invariants_paired_end_matrix() {
    let dir = tempdir().unwrap();
    let fq1 = dir.path().join("pe_r1.fastq");
    let fq2 = dir.path().join("pe_r2.fastq");
    let bc = dir.path().join("barcodes.csv");

    // Dual barcode PE:
    // sample1: BC1=AAGTCCAA, BC2=GGAGTACT
    // sample2: BC1=TTTTTTTT, BC2=CCCCCCCC
    fs::write(
        &bc,
        "SampleNumber,Barcode1,Barcode2\n\
         sample1,AAGTCCAA,GGAGTACT\n\
         sample2,TTTTTTTT,CCCCCCCC\n",
    )
    .unwrap();

    // 8 read pairs:
    // 1. Long canonical match to sample1 (R1=AAGTCCAA..., R2=GGAGTACT...)
    // 2. Short canonical match to sample1 (R1 trimmed < 10)
    // 3. Long swapped match to sample1 (R1=GGAGTACT..., R2=AAGTCCAA...)
    // 4. Short swapped match to sample1 (R2 trimmed < 10)
    // 5. Long ambiguous pair
    // 6. Short ambiguous pair
    // 7. Long no-match pair
    // 8. Short no-match pair
    write_fastq(
        &fq1,
        &[
            (
                "p1_long_canon",
                "AAGTCCAACCCCCCCCCCCC",
                "IIIIIIIIIIIIIIIIIIII",
            ),
            ("p2_short_canon", "AAGTCCAAACCC", "IIIIIIIIIIII"),
            (
                "p3_long_swap",
                "GGAGTACTCCCCCCCCCCCC",
                "IIIIIIIIIIIIIIIIIIII",
            ),
            (
                "p4_short_swap",
                "GGAGTACTCCCCCCCCCCCC",
                "IIIIIIIIIIIIIIIIIIII",
            ),
            (
                "p5_long_ambig",
                "AAGTTTTTCCCCCCCCCCCC",
                "IIIIIIIIIIIIIIIIIIII",
            ),
            ("p6_short_ambig", "AAGTTTTTACCC", "IIIIIIIIIIII"),
            (
                "p7_long_nomatch",
                "NNNNNNNNCCCCCCCCCCCC",
                "IIIIIIIIIIIIIIIIIIII",
            ),
            ("p8_short_nomatch", "NNNNNNNNACCC", "IIIIIIIIIIII"),
        ],
    );
    write_fastq(
        &fq2,
        &[
            (
                "p1_long_canon",
                "GGAGTACTCCCCCCCCCCCC",
                "IIIIIIIIIIIIIIIIIIII",
            ),
            (
                "p2_short_canon",
                "GGAGTACTCCCCCCCCCCCC",
                "IIIIIIIIIIIIIIIIIIII",
            ),
            (
                "p3_long_swap",
                "AAGTCCAACCCCCCCCCCCC",
                "IIIIIIIIIIIIIIIIIIII",
            ),
            ("p4_short_swap", "AAGTCCAAACCC", "IIIIIIIIIIII"),
            (
                "p5_long_ambig",
                "GGAGCCCCCCCCCCCCCCCC",
                "IIIIIIIIIIIIIIIIIIII",
            ),
            ("p6_short_ambig", "GGAGCCCCACCC", "IIIIIIIIIIII"),
            (
                "p7_long_nomatch",
                "NNNNNNNNCCCCCCCCCCCC",
                "IIIIIIIIIIIIIIIIIIII",
            ),
            ("p8_short_nomatch", "NNNNNNNNACCC", "IIIIIIIIIIII"),
        ],
    );

    for min_length in [0, 10] {
        for threads in [1, 2] {
            let out_dir = dir.path().join(format!("out_pe_l{min_length}_t{threads}"));

            Command::cargo_bin("seqmux")
                .unwrap()
                .args([
                    "demux",
                    "-i",
                    fq1.to_str().unwrap(),
                    "-I",
                    fq2.to_str().unwrap(),
                    "-b",
                    bc.to_str().unwrap(),
                    "-o",
                    out_dir.to_str().unwrap(),
                    "--min-length",
                    &min_length.to_string(),
                    "--threads",
                    &threads.to_string(),
                    "--mismatches-1",
                    "4",
                    "--mismatches-2",
                    "4",
                    "--no-adapter",
                    "--no-gzip",
                ])
                .assert()
                .success();

            let summary_file = out_dir.join("seqmux.summary.tsv");
            let content = fs::read_to_string(&summary_file).unwrap();
            let (lines, map) = parse_summary_tsv(&content);

            check_invariants(&map, true);

            let total: u64 = map["total_reads"].parse().unwrap();
            assert_eq!(total, 8);

            let matched_before: u64 = map["matched_before_filter"].parse().unwrap();
            let canon: u64 = map["orientation_canonical"].parse().unwrap();
            let swap: u64 = map["orientation_swapped"].parse().unwrap();
            assert_eq!(matched_before, 4);
            assert_eq!(canon, 2);
            assert_eq!(swap, 2);
            assert_eq!(canon + swap, matched_before);

            // TSV suffix check
            let n = lines.len();
            assert!(lines[n - 4].starts_with("matched_before_filter\t"));
            assert!(lines[n - 3].starts_with("no_match_before_filter\t"));
            assert!(lines[n - 2].starts_with("ambiguous_after_filter\t"));
            assert!(lines[n - 1].starts_with("assignment_rate_before_filter\t"));
        }
    }
}
