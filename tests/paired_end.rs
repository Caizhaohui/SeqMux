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
