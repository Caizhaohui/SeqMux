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

fn setup_fixtures(
    dir: &std::path::Path,
) -> (std::path::PathBuf, std::path::PathBuf, std::path::PathBuf) {
    let fq1 = dir.join("r1.fastq");
    let fq2 = dir.join("r2.fastq");
    let bc = dir.join("barcodes.csv");

    write_fastq(&fq1, &[("read1", "AAGTCCAACCCCCCCC", "IIIIIIIIIIIIIIII")]);
    write_fastq(&fq2, &[("read1", "GGAGTACTCCCCCCCC", "IIIIIIIIIIIIIIII")]);
    fs::write(
        &bc,
        "SampleNumber,Barcode1,Barcode2\n\
         sample1,AAGTCCAA,GGAGTACT\n",
    )
    .unwrap();

    (fq1, fq2, bc)
}

#[test]
fn test_boundary_threads_zero_rejected() {
    let dir = tempdir().unwrap();
    let (fq1, _fq2, bc) = setup_fixtures(dir.path());
    let out = dir.path().join("out_nonexistent");

    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--threads",
            "0",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--threads must be at least 1"));

    assert!(
        !out.exists(),
        "out directory must not be created on validation failure"
    );
}

#[test]
fn test_boundary_threads_one_accepted() {
    let dir = tempdir().unwrap();
    let (fq1, _fq2, bc) = setup_fixtures(dir.path());
    let out = dir.path().join("out");

    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--threads",
            "1",
        ])
        .assert()
        .success();

    assert!(out.exists());
}

#[test]
fn test_boundary_chunk_reads_zero_rejected() {
    let dir = tempdir().unwrap();
    let (fq1, _fq2, bc) = setup_fixtures(dir.path());
    let out = dir.path().join("out_nonexistent");

    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--chunk-reads",
            "0",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--chunk-reads must be at least 1"));

    assert!(
        !out.exists(),
        "out directory must not be created on validation failure"
    );
}

#[test]
fn test_boundary_chunk_reads_one_accepted() {
    let dir = tempdir().unwrap();
    let (fq1, _fq2, bc) = setup_fixtures(dir.path());
    let out = dir.path().join("out");

    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--chunk-reads",
            "1",
        ])
        .assert()
        .success();

    assert!(out.exists());
}

#[test]
fn test_boundary_compression_level_boundaries() {
    let dir = tempdir().unwrap();
    let (fq1, _fq2, bc) = setup_fixtures(dir.path());

    // level 0 rejected
    let out0 = dir.path().join("out0");
    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out0.to_str().unwrap(),
            "--compression-level",
            "0",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "--compression-level must be between 1 and 9",
        ));
    assert!(!out0.exists());

    // level 10 rejected
    let out10 = dir.path().join("out10");
    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out10.to_str().unwrap(),
            "--compression-level",
            "10",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "--compression-level must be between 1 and 9",
        ));
    assert!(!out10.exists());

    // level 1 accepted
    let out1 = dir.path().join("out1");
    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out1.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .assert()
        .success();
    assert!(out1.exists());

    // level 9 accepted
    let out9 = dir.path().join("out9");
    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out9.to_str().unwrap(),
            "--compression-level",
            "9",
        ])
        .assert()
        .success();
    assert!(out9.exists());
}

#[test]
fn test_boundary_adapter_error_rate_boundaries() {
    let dir = tempdir().unwrap();
    let (fq1, _fq2, bc) = setup_fixtures(dir.path());

    // negative rejected
    let out_neg = dir.path().join("out_neg");
    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out_neg.to_str().unwrap(),
            "--adapter-error-rate",
            "-0.1",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--adapter-error-rate"));
    assert!(!out_neg.exists());

    // > 1.0 rejected
    let out_high = dir.path().join("out_high");
    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out_high.to_str().unwrap(),
            "--adapter-error-rate",
            "1.01",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--adapter-error-rate"));
    assert!(!out_high.exists());

    // NaN rejected
    let out_nan = dir.path().join("out_nan");
    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out_nan.to_str().unwrap(),
            "--adapter-error-rate",
            "NaN",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--adapter-error-rate"));
    assert!(!out_nan.exists());

    // inf rejected
    let out_inf = dir.path().join("out_inf");
    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out_inf.to_str().unwrap(),
            "--adapter-error-rate",
            "inf",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--adapter-error-rate"));
    assert!(!out_inf.exists());

    // 0.0 accepted
    let out0 = dir.path().join("out_aer0");
    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out0.to_str().unwrap(),
            "--adapter-error-rate",
            "0.0",
        ])
        .assert()
        .success();
    assert!(out0.exists());

    // 1.0 accepted
    let out1 = dir.path().join("out_aer1");
    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out1.to_str().unwrap(),
            "--adapter-error-rate",
            "1.0",
        ])
        .assert()
        .success();
    assert!(out1.exists());
}

#[test]
fn test_boundary_min_adapter_overlap_boundaries() {
    let dir = tempdir().unwrap();
    let (fq1, _fq2, bc) = setup_fixtures(dir.path());

    // 0 rejected
    let out0 = dir.path().join("out_mao0");
    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out0.to_str().unwrap(),
            "--min-adapter-overlap",
            "0",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "--min-adapter-overlap must be at least 1",
        ));
    assert!(!out0.exists());

    // 1 accepted
    let out1 = dir.path().join("out_mao1");
    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out1.to_str().unwrap(),
            "--min-adapter-overlap",
            "1",
        ])
        .assert()
        .success();
    assert!(out1.exists());
}

#[test]
fn test_input_and_input2_same_file_rejected() {
    let dir = tempdir().unwrap();
    let (fq1, _fq2, bc) = setup_fixtures(dir.path());
    let out = dir.path().join("out_same");

    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq1.to_str().unwrap(),
            "-I",
            fq1.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "input and input2 cannot point to the same file",
        ));

    assert!(!out.exists(), "out directory must not be created");
}

#[test]
fn test_validate_subcommand_input_and_input2_same_file_rejected() {
    let dir = tempdir().unwrap();
    let (fq1, _fq2, _bc) = setup_fixtures(dir.path());

    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "validate",
            "-i",
            fq1.to_str().unwrap(),
            "-I",
            fq1.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "input and input2 cannot point to the same file",
        ));
}
