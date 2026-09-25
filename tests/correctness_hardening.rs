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

// -------------------------------------------------------------------------
// 1. Sanitized sample-name output collisions
// -------------------------------------------------------------------------

#[test]
fn test_cli_rejects_sanitized_sample_name_collision() {
    let dir = tempdir().unwrap();
    let fq = dir.path().join("r.fastq");
    let bc = dir.path().join("barcodes.csv");
    let out = dir.path().join("out");

    write_fastq(&fq, &[("r1", "AAGTCCAACCCCCCCC", "IIIIIIIIIIIIIIII")]);
    // "sample:1" and "sample/1" both sanitize to "sample_1"
    fs::write(
        &bc,
        "SampleNumber,Barcode1,Barcode2\n\
         sample:1,AAGTCCAA,GGAGTACT\n\
         sample/1,GACCTGAA,GGACTTGG\n",
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
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("collides"));
}

#[test]
fn test_cli_rejects_unassigned_sample_name() {
    let dir = tempdir().unwrap();
    let fq = dir.path().join("r.fastq");
    let bc = dir.path().join("barcodes.csv");
    let out = dir.path().join("out");

    write_fastq(&fq, &[("r1", "AAGTCCAACCCCCCCC", "IIIIIIIIIIIIIIII")]);
    fs::write(
        &bc,
        "SampleNumber,Barcode1,Barcode2\n\
         unassigned,AAGTCCAA,GGAGTACT\n",
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
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "reserved output name 'unassigned'",
        ));
}

// -------------------------------------------------------------------------
// 2. Prefix/path collisions
// -------------------------------------------------------------------------

#[test]
fn test_cli_rejects_empty_or_whitespace_prefix() {
    let dir = tempdir().unwrap();
    let fq = dir.path().join("r.fastq");
    let bc = dir.path().join("barcodes.csv");

    write_fastq(&fq, &[("r1", "AAGTCCAACCCCCCCC", "IIIIIIIIIIIIIIII")]);
    fs::write(&bc, "SampleNumber,Barcode1\ns1,AAGTCCAA\n").unwrap();

    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-p",
            "",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("output prefix cannot be empty"));
}

#[test]
fn test_cli_rejects_path_separator_in_prefix() {
    let dir = tempdir().unwrap();
    let fq = dir.path().join("r.fastq");
    let bc = dir.path().join("barcodes.csv");

    write_fastq(&fq, &[("r1", "AAGTCCAACCCCCCCC", "IIIIIIIIIIIIIIII")]);
    fs::write(&bc, "SampleNumber,Barcode1\ns1,AAGTCCAA\n").unwrap();

    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-p",
            "sub/dir",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot contain path separators"));
}

#[test]
fn test_cli_rejects_identical_inputs() {
    let dir = tempdir().unwrap();
    let fq = dir.path().join("r.fastq");
    let bc = dir.path().join("barcodes.csv");

    write_fastq(&fq, &[("r1", "AAGTCCAACCCCCCCC", "IIIIIIIIIIIIIIII")]);
    fs::write(
        &bc,
        "SampleNumber,Barcode1,Barcode2\ns1,AAGTCCAA,GGAGTACT\n",
    )
    .unwrap();

    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            fq.to_str().unwrap(),
            "-I",
            fq.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot point to the same file"));
}

#[test]
fn test_cli_rejects_output_overwriting_input() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("out");
    fs::create_dir_all(&out).unwrap();
    let fq = out.join("test_s1.fastq");
    let bc = dir.path().join("barcodes.csv");

    write_fastq(&fq, &[("r1", "AAGTCCAACCCCCCCC", "IIIIIIIIIIIIIIII")]);
    fs::write(&bc, "SampleNumber,Barcode1\ns1,AAGTCCAA\n").unwrap();

    // Planned output is out/test_s1.fastq which is the input file!
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
            "--no-gzip",
            "--force", // Even with force, input overwrite is fatal
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("would overwrite input file"));

    // Ensure input file was not wiped/truncated
    let content = fs::read_to_string(&fq).unwrap();
    assert!(content.contains("@r1"));
}

// -------------------------------------------------------------------------
// 3. Output preflight and --force semantics
// -------------------------------------------------------------------------

#[test]
fn test_preflight_blocks_before_reading_when_output_exists_without_force() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("out");
    fs::create_dir_all(&out).unwrap();

    let fq = dir.path().join("r.fastq");
    let bc = dir.path().join("barcodes.csv");

    write_fastq(&fq, &[("r1", "AAGTCCAACCCCCCCC", "IIIIIIIIIIIIIIII")]);
    fs::write(
        &bc,
        "SampleNumber,Barcode1\n\
         s1,AAGTCCAA\n\
         s2,GACCTGAA\n",
    )
    .unwrap();

    // Pre-create output for s2 (which won't even be encountered in the reads!)
    let existing_s2 = out.join("test_s2.fastq");
    fs::write(&existing_s2, b"pre-existing s2").unwrap();

    // Run without --force
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
            "--no-gzip",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("already exists"));

    // Verify s1 output was NOT created (preflight stopped execution before writing anything)
    assert!(!out.join("test_s1.fastq").exists());
    // Verify existing_s2 was NOT overwritten
    assert_eq!(fs::read_to_string(&existing_s2).unwrap(), "pre-existing s2");

    // Now run WITH --force: should succeed and overwrite
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
            "--no-gzip",
            "--force",
        ])
        .assert()
        .success();

    assert!(out.join("test_s1.fastq").exists());
}

// -------------------------------------------------------------------------
// 4. Summary overwrite semantics
// -------------------------------------------------------------------------

#[test]
fn test_summary_overwrite_fails_without_force() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("out");
    fs::create_dir_all(&out).unwrap();

    let fq = dir.path().join("r.fastq");
    let bc = dir.path().join("barcodes.csv");
    let summary = out.join("my_summary.tsv");

    write_fastq(&fq, &[("r1", "AAGTCCAACCCCCCCC", "IIIIIIIIIIIIIIII")]);
    fs::write(&bc, "SampleNumber,Barcode1\ns1,AAGTCCAA\n").unwrap();
    fs::write(&summary, "existing summary content").unwrap();

    // Normal mode without force
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
            "--summary",
            summary.to_str().unwrap(),
            "--no-gzip",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("already exists"));

    assert_eq!(
        fs::read_to_string(&summary).unwrap(),
        "existing summary content"
    );

    // Counts-only mode without force
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
            "--summary",
            summary.to_str().unwrap(),
            "--counts-only",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("already exists"));

    assert_eq!(
        fs::read_to_string(&summary).unwrap(),
        "existing summary content"
    );

    // Counts-only mode WITH force
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
            "--summary",
            summary.to_str().unwrap(),
            "--counts-only",
            "--force",
        ])
        .assert()
        .success();

    let new_content = fs::read_to_string(&summary).unwrap();
    assert!(new_content.contains("total_reads\t1"));
}

// -------------------------------------------------------------------------
// 6. PE pair-vs-mate statistics semantics
// -------------------------------------------------------------------------

#[test]
fn test_pe_statistics_count_pairs_not_mates() {
    let dir = tempdir().unwrap();
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let bc = dir.path().join("barcodes.csv");
    let out = dir.path().join("out");

    // Illumina adapters:
    // R1: AGATCGGAAGAGCACACGTCTGAA
    // R2: AGATCGGAAGAGCGTCGTG
    // 1 pair: both R1 and R2 have 3' adapter AND low quality bases at 3'
    write_fastq(
        &r1,
        &[(
            "read1/1",
            "AAGTCCAAAAAAAAAAAGATCGGAAGAGCACACGTCTGAA",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!!!!!!",
        )],
    );
    write_fastq(
        &r2,
        &[(
            "read1/2",
            "GGAGTACTGGGGGGGGAGATCGGAAGAGCGTCGTG",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII!!!!!",
        )],
    );

    fs::write(
        &bc,
        "SampleNumber,Barcode1,Barcode2\n\
         I464469-A1,AAGTCCAA,GGAGTACT\n",
    )
    .unwrap();

    Command::cargo_bin("seqmux")
        .unwrap()
        .args([
            "demux",
            "-i",
            r1.to_str().unwrap(),
            "-I",
            r2.to_str().unwrap(),
            "-b",
            bc.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "-p",
            "pe_stat",
            "-q",
            "20",
            "--no-gzip",
            "--force",
        ])
        .assert()
        .success();

    let summary = fs::read_to_string(out.join("pe_stat.summary.tsv")).unwrap();

    // Total reads in PE is 1 pair
    assert!(summary.contains("total_reads\t1"), "summary: {summary}");
    // Both R1 and R2 had adapter trimmed, but adapter_trimmed must be 1 (pairs), not 2
    assert!(
        summary.contains("adapter_trimmed\t1"),
        "summary should report 1 adapter_trimmed pair: {summary}"
    );
    // Both R1 and R2 had quality trimmed, but quality_trimmed must be 1 (pairs), not 2
    assert!(
        summary.contains("quality_trimmed\t1"),
        "summary should report 1 quality_trimmed pair: {summary}"
    );
}
