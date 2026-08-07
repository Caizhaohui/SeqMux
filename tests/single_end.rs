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

    // Dual barcode SE: Barcode1 at 5', Barcode2 at 3'
    write_fastq(
        &fq,
        &[
            ("r1", "AAGTCCAAAAAAAAAAGGAGTACT", "IIIIIIIIIIIIIIIIIIIIIIII"),
            ("r2", "GACCTGAAAAAAAAAAGGACTTGG", "IIIIIIIIIIIIIIIIIIIIIIII"),
            ("r3", "GGGGGGGGGGGGGGGGGGGGGGGG", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ],
    );
    fs::write(
        &bc,
        "SampleNumber,Barcode1,Barcode2\n\
         sampleA,AAGTCCAA,GGAGTACT\n\
         sampleB,GACCTGAA,GGACTTGG\n",
    )
    .unwrap();

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
    // both barcodes trimmed: 24 - 8 - 8 = 8
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

    // Barcode1 ATGATGAT with 1 mismatch → ATGATGA A
    write_fastq(&fq, &[("r1", "ATGATGAAAAAAAAA", "IIIIIIIIIIIIIII")]);
    fs::write(
        &bc,
        "SampleNumber,Barcode1\n\
         s1,ATGATGAT\n\
         s2,CCGTAACG\n",
    )
    .unwrap();

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
            "--mismatches-1",
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

#[test]
fn reject_ultraplex_format() {
    let dir = tempdir().unwrap();
    let fq = dir.path().join("reads.fastq");
    let bc = dir.path().join("barcodes.csv");
    let out = dir.path().join("out");
    write_fastq(&fq, &[("r1", "ATGAAAAA", "IIIIIIII")]);
    // Old Ultraplex format without proper header
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
            "x",
            "-t",
            "1",
            "--no-gzip",
            "--force",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("SampleNumber").or(predicate::str::contains("header")));
}
