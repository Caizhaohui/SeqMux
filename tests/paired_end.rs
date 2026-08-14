use assert_cmd::Command;
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
fn paired_dual_barcode_demux() {
    let dir = tempdir().unwrap();
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let bc = dir.path().join("barcodes.csv");
    let out = dir.path().join("out");

    // Barcode1 on R1 5', Barcode2 on R2 5'
    write_fastq(
        &r1,
        &[
            ("readA/1", "AAGTCCAACCCCCCCC", "IIIIIIIIIIIIIIII"),
            ("readB/1", "GGGGGGGGCCCCCCCC", "IIIIIIIIIIIIIIII"),
        ],
    );
    write_fastq(
        &r2,
        &[
            ("readA/2", "GGAGTACTGGGGGGGG", "IIIIIIIIIIIIIIII"),
            ("readB/2", "GGGGGGGGGGGGGGGG", "IIIIIIIIIIIIIIII"),
        ],
    );
    fs::write(
        &bc,
        "SampleNumber,Barcode1,Barcode2,PCR_product,rawdata1,rawdata2,library_round,F_primer,R_primer\n\
         I464469-A1,AAGTCCAA,GGAGTACT,atgc,x_1.fq.gz,x_2.fq.gz,A,fwd,rev\n",
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
            "pe",
            "-t",
            "2",
            "--no-gzip",
            "--force",
        ])
        .assert()
        .success();

    let o1 = fs::read_to_string(out.join("pe_I464469-A1_R1.fastq")).unwrap();
    assert!(o1.contains("@readA"));
    assert!(o1.contains("\nCCCCCCCC\n")); // Barcode1 trimmed

    let o2 = fs::read_to_string(out.join("pe_I464469-A1_R2.fastq")).unwrap();
    assert!(o2.contains("@readA"));
    assert!(o2.contains("\nGGGGGGGG\n")); // Barcode2 trimmed

    let u1 = fs::read_to_string(out.join("pe_unassigned_R1.fastq")).unwrap();
    assert!(u1.contains("@readB"));
}

#[test]
fn paired_swapped_orientation_is_assigned_and_canonicalized() {
    let dir = tempdir().unwrap();
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let bc = dir.path().join("barcodes.csv");
    let out = dir.path().join("out");

    // R1 has Barcode2, R2 has Barcode1 (swapped vs sample table)
    write_fastq(&r1, &[("swapA/1", "GGAGTACTGGGGGGGG", "IIIIIIIIIIIIIIII")]);
    write_fastq(&r2, &[("swapA/2", "AAGTCCAACCCCCCCC", "IIIIIIIIIIIIIIII")]);
    fs::write(
        &bc,
        "SampleNumber,Barcode1,Barcode2\nI464469-A1,AAGTCCAA,GGAGTACT\n",
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
            "sw",
            "-t",
            "1",
            "--no-gzip",
            "--force",
        ])
        .assert()
        .success();

    let o1 = fs::read_to_string(out.join("sw_I464469-A1_R1.fastq")).unwrap();
    // Canonicalize: output R1 is the mate that had Barcode1, then trimmed
    assert!(o1.contains("@swapA/2"));
    assert!(o1.contains("\nCCCCCCCC\n"));

    let o2 = fs::read_to_string(out.join("sw_I464469-A1_R2.fastq")).unwrap();
    assert!(o2.contains("@swapA/1"));
    assert!(o2.contains("\nGGGGGGGG\n"));
}

#[test]
fn paired_orientation_canonical_drops_swapped() {
    let dir = tempdir().unwrap();
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let bc = dir.path().join("barcodes.csv");
    let out = dir.path().join("out");

    write_fastq(&r1, &[("swapA/1", "GGAGTACTGGGGGGGG", "IIIIIIIIIIIIIIII")]);
    write_fastq(&r2, &[("swapA/2", "AAGTCCAACCCCCCCC", "IIIIIIIIIIIIIIII")]);
    fs::write(
        &bc,
        "SampleNumber,Barcode1,Barcode2\nI464469-A1,AAGTCCAA,GGAGTACT\n",
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
            "can",
            "-t",
            "1",
            "--no-gzip",
            "--orientation",
            "canonical",
            "--force",
        ])
        .assert()
        .success();

    let u1 = fs::read_to_string(out.join("can_unassigned_R1.fastq")).unwrap();
    assert!(u1.contains("@swapA"));
    assert!(!out.join("can_I464469-A1_R1.fastq").exists());
}

#[test]
fn i464_real_fixture_both_orientations() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("out");
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let r1 = root.join("i464_real_R1.fastq");
    let r2 = root.join("i464_real_R2.fastq");
    let bc = root.join("I464-469erdai_barcode_and_name.csv");

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
            "i464",
            "-t",
            "2",
            "--no-gzip",
            "--force",
        ])
        .assert()
        .success();

    let b2 = fs::read_to_string(out.join("i464_I464469-B2_R1.fastq")).unwrap();
    assert!(b2.contains("00013998"));

    let b6 = fs::read_to_string(out.join("i464_I464469-B6_R1.fastq")).unwrap();
    assert!(b6.contains("00014336"));

    // swapped C5: after canonicalize, R1 should be the Barcode1/F-primer mate
    let c5_r1 = fs::read_to_string(out.join("i464_I464469-C5_R1.fastq")).unwrap();
    assert!(c5_r1.contains("00010916"));
    assert!(c5_r1.contains("\nACCAAATGCATTCGCATTGCG"));

    let a4 = fs::read_to_string(out.join("i464_I464469-A4_R1.fastq")).unwrap();
    assert!(a4.contains("00014624"));

    let u = fs::read_to_string(out.join("i464_unassigned_R1.fastq")).unwrap();
    assert!(u.contains("00011598"));

    let summary = fs::read_to_string(out.join("i464.summary.tsv")).unwrap();
    assert!(summary.contains("orientation_canonical\t2"));
    assert!(summary.contains("orientation_swapped\t2"));
    assert!(summary.contains("assigned\t4"));
    assert!(summary.contains("unassigned\t1"));
}

#[test]
fn paired_id_mismatch_errors() {
    let dir = tempdir().unwrap();
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let bc = dir.path().join("barcodes.csv");
    let out = dir.path().join("out");

    write_fastq(&r1, &[("readA/1", "AAGTCCAACCCC", "IIIIIIIIIIII")]);
    write_fastq(&r2, &[("readZ/2", "GGAGTACTGGGG", "IIIIIIIIIIII")]);
    fs::write(
        &bc,
        "SampleNumber,Barcode1,Barcode2\nI464469-A1,AAGTCCAA,GGAGTACT\n",
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
            "bad",
            "-t",
            "1",
            "--no-gzip",
            "--force",
        ])
        .assert()
        .failure();
}
