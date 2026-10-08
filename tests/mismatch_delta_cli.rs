use assert_cmd::Command;
use predicates::prelude::*;
use std::fs;
use tempfile::tempdir;

#[test]
fn margin_rejections_use_existing_unassigned_and_summary_paths() {
    let dir = tempdir().unwrap();
    let reads = dir.path().join("reads.fastq");
    let barcodes = dir.path().join("barcodes.csv");
    fs::write(&reads, "@read1\nAAAATTTT\n+\nIIIIIIII\n").unwrap();
    fs::write(&barcodes, "SampleNumber,Barcode1\nA,AAAC\nB,AACC\n").unwrap();

    for (label, flags, assigned, ambiguous) in [
        ("default", vec![], 1, 0),
        ("zero", vec!["--min-mismatch-delta", "0"], 1, 0),
        ("margin", vec!["--min-mismatch-delta", "2"], 0, 1),
        (
            "discard",
            vec!["--min-mismatch-delta", "2", "--discard-unassigned"],
            0,
            1,
        ),
        (
            "counts",
            vec!["--min-mismatch-delta", "2", "--counts-only"],
            0,
            1,
        ),
    ] {
        let out = dir.path().join(label);
        Command::cargo_bin("seqmux")
            .unwrap()
            .arg("demux")
            .arg("-i")
            .arg(&reads)
            .arg("-b")
            .arg(&barcodes)
            .arg("-o")
            .arg(&out)
            .args([
                "--mismatches-1",
                "2",
                "--no-adapter",
                "--no-gzip",
                "-t",
                "2",
            ])
            .args(flags)
            .assert()
            .success();

        let summary = fs::read_to_string(out.join("seqmux.summary.tsv")).unwrap();
        for line in [
            "total_reads\t1".to_string(),
            format!("assigned\t{assigned}"),
            format!("unassigned\t{ambiguous}"),
            format!("ambiguous\t{ambiguous}"),
        ] {
            assert!(summary.lines().any(|actual| actual == line), "{summary}");
        }
        if assigned == 1 {
            assert_eq!(
                fs::read_to_string(out.join("seqmux_A.fastq")).unwrap(),
                "@read1\nTTTT\n+\nIIII\n"
            );
        } else {
            assert!(!out.join("seqmux_A.fastq").exists());
        }
        if label == "margin" {
            assert_eq!(
                fs::read_to_string(out.join("seqmux_unassigned.fastq")).unwrap(),
                fs::read_to_string(&reads).unwrap()
            );
        } else {
            assert!(!out.join("seqmux_unassigned.fastq").exists());
        }
    }
}

#[test]
fn paired_low_margin_keeps_mates_in_unassigned_output() {
    let dir = tempdir().unwrap();
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let barcodes = dir.path().join("barcodes.csv");
    let out = dir.path().join("out");
    fs::write(&r1, "@pair/1\nAAAATTTT\n+\nIIIIIIII\n").unwrap();
    fs::write(&r2, "@pair/2\nCCCCGGGG\n+\nIIIIIIII\n").unwrap();
    fs::write(
        &barcodes,
        "SampleNumber,Barcode1,Barcode2\nA,AAAC,CCCC\nB,AACA,CCCA\n",
    )
    .unwrap();
    Command::cargo_bin("seqmux")
        .unwrap()
        .arg("demux")
        .arg("-i")
        .arg(&r1)
        .arg("-I")
        .arg(&r2)
        .arg("-b")
        .arg(&barcodes)
        .arg("-o")
        .arg(&out)
        .args([
            "--mismatches-1",
            "1",
            "--mismatches-2",
            "1",
            "--min-mismatch-delta",
            "2",
            "--no-adapter",
            "--no-gzip",
            "-t",
            "2",
        ])
        .assert()
        .success();
    let summary = fs::read_to_string(out.join("seqmux.summary.tsv")).unwrap();
    for line in [
        "total_reads\t1",
        "assigned\t0",
        "unassigned\t1",
        "ambiguous\t1",
        "orientation_canonical\t0",
        "orientation_swapped\t0",
    ] {
        assert!(summary.lines().any(|actual| actual == line), "{summary}");
    }
    for (mate, input) in [("R1", r1), ("R2", r2)] {
        assert_eq!(
            fs::read(out.join(format!("seqmux_unassigned_{mate}.fastq"))).unwrap(),
            fs::read(input).unwrap()
        );
    }
}

#[test]
fn positive_margin_preserves_exact_i464_fixture_outputs_and_summary() {
    let dir = tempdir().unwrap();
    let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    for (label, flags) in [
        ("default", vec![]),
        ("positive", vec!["--min-mismatch-delta", "2"]),
    ] {
        Command::cargo_bin("seqmux")
            .unwrap()
            .arg("demux")
            .arg("-i")
            .arg(fixtures.join("i464_real_R1.fastq"))
            .arg("-I")
            .arg(fixtures.join("i464_real_R2.fastq"))
            .arg("-b")
            .arg(fixtures.join("I464-469erdai_barcode_and_name.csv"))
            .arg("-o")
            .arg(dir.path().join(label))
            .args(["--no-gzip", "-t", "2"])
            .args(flags)
            .assert()
            .success();
    }
    let filenames = |label| {
        let mut names: Vec<_> = fs::read_dir(dir.path().join(label))
            .unwrap()
            .map(|entry| entry.unwrap().file_name())
            .collect();
        names.sort();
        names
    };
    assert_eq!(filenames("default"), filenames("positive"));
    for name in filenames("default") {
        assert_eq!(
            fs::read(dir.path().join("default").join(&name)).unwrap(),
            fs::read(dir.path().join("positive").join(&name)).unwrap(),
            "{name:?}"
        );
    }
}

#[test]
fn margin_cli_help_and_invalid_values() {
    Command::cargo_bin("seqmux")
        .unwrap()
        .args(["demux", "--help"])
        .assert()
        .success()
        .stdout(predicate::str::contains("--min-mismatch-delta <INT>"))
        .stdout(predicate::str::contains(
            "closest distinct competitor globally",
        ));
    for value in ["-1", "1.5", "18446744073709551616"] {
        Command::cargo_bin("seqmux")
            .unwrap()
            .args(["demux", "-i", "unused.fastq", "-b", "unused.csv"])
            .arg(format!("--min-mismatch-delta={value}"))
            .assert()
            .failure()
            .stderr(predicate::str::contains("--min-mismatch-delta"));
    }
}
