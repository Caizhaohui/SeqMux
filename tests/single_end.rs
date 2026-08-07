use assert_cmd::Command;
use predicates::prelude::*;
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

#[test]
fn demux_two_samples_exact() {
    let dir = tempdir().unwrap();
    let fq = dir.path().join("reads.fastq");
    let bc = dir.path().join("barcodes.csv");
    let out = dir.path().join("out");

    // NNNATGNN:sampleA  NNNCCGNN:sampleB
    write_fastq(
        &fq,
        &[
            ("r1", "ACGATGTCAAAAAAAA", "IIIIIIIIIIIIIIII"), // UMI ACGTC, bc ATG
            ("r2", "TTTCCGNNAAAAAAAA", "IIIIIIIIIIIIIIII"), // UMI TTTNA? wait NNNCCGNN -> TTTCCGNN
            ("r3", "GGGAAAGGGAAAAAAA", "IIIIIIIIIIIIIIII"), // no match
        ],
    );
    // Fix r2: NNNCCGNN positions: 0-2 UMI, 3-5 CCG, 6-7 UMI
    // seq TTTCCGTAAAAAAA → UMI TTTAA, bc CCG
    write_fastq(
        &fq,
        &[
            ("r1", "ACGATGTCAAAAAAAA", "IIIIIIIIIIIIIIII"),
            ("r2", "TTTCCGTAAAAAAA", "IIIIIIIIIIIIII"),
            ("r3", "GGGAAAGGGAAAAAAA", "IIIIIIIIIIIIIIII"),
        ],
    );
    fs::write(&bc, "NNNATGNN:sampleA\nNNNCCGNN:sampleB\n").unwrap();

    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "-p",
            "test",
            "-t",
            "1",
            "--no-gzip",
            "--force",
        ])
        .assert()
        .success();

    let a = fs::read_to_string(out.join("test_sampleA.fastq")).unwrap();
    assert!(a.contains("@r1"));
    assert!(a.contains("rbc:ACGTC"));
    // barcode trimmed: original ACGATGTCAAAAAAAA len 16, trim 8 → AAAAAAAA
    assert!(a.contains("\nAAAAAAAA\n"));

    let b = fs::read_to_string(out.join("test_sampleB.fastq")).unwrap();
    assert!(b.contains("@r2"));

    let u = fs::read_to_string(out.join("test_unassigned.fastq")).unwrap();
    assert!(u.contains("@r3"));
}

#[test]
fn validate_command() {
    let dir = tempdir().unwrap();
    let fq = dir.path().join("reads.fastq");
    write_fastq(&fq, &[("r1", "ACGT", "IIII")]);

    Command::cargo_bin("seqmux")
        .unwrap()
        .args(["validate", "-i", fq.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains("records\t1"));
}

#[test]
fn mismatch_tolerance() {
    let dir = tempdir().unwrap();
    let fq = dir.path().join("reads.fastq");
    let bc = dir.path().join("barcodes.csv");
    let out = dir.path().join("out");

    // ATG with 1 mismatch → AAG
    write_fastq(&fq, &[("r1", "AAGAAAAA", "IIIIIIII")]);
    fs::write(&bc, "ATG:s1\nCCG:s2\n").unwrap();

    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "-p",
            "m",
            "-t",
            "1",
            "--no-gzip",
            "--mismatches-5",
            "1",
            "--force",
        ])
        .assert()
        .success();

    let a = fs::read_to_string(out.join("m_s1.fastq")).unwrap();
    assert!(a.contains("@r1"));
}

#[test]
fn help_and_version() {
    Command::cargo_bin("seqmux")
        .unwrap()
        .arg("--help")
        .assert()
        .success()
        .stdout(predicate::str::contains("SeqMux"));

    Command::cargo_bin("seqmux")
        .unwrap()
        .arg("--version")
        .assert()
        .success()
        .stdout(predicate::str::contains("seqmux"));
}
