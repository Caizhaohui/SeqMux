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
fn paired_demux_basic() {
    let dir = tempdir().unwrap();
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let bc = dir.path().join("barcodes.csv");
    let out = dir.path().join("out");

    // 5' barcode ATG on R1
    write_fastq(
        &r1,
        &[
            ("readA/1", "ATGCCCCCCCC", "IIIIIIIIIII"),
            ("readB/1", "GGGCCCCCCCC", "IIIIIIIIIII"),
        ],
    );
    write_fastq(
        &r2,
        &[
            ("readA/2", "GGGGGGCAT", "IIIIIIIII"), // RC(ATG)=CAT at end
            ("readB/2", "GGGGGGGGG", "IIIIIIIII"),
        ],
    );
    fs::write(&bc, "ATG:sampleX\n").unwrap();

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

    let o1 = fs::read_to_string(out.join("pe_sampleX_R1.fastq")).unwrap();
    assert!(o1.contains("@readA"));
    assert!(o1.contains("\nCCCCCCCC\n")); // ATG trimmed

    let o2 = fs::read_to_string(out.join("pe_sampleX_R2.fastq")).unwrap();
    assert!(o2.contains("@readA"));

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

    write_fastq(&r1, &[("readA/1", "ATGCCC", "IIIIII")]);
    write_fastq(&r2, &[("readZ/2", "GGGCAT", "IIIIII")]);
    fs::write(&bc, "ATG:s1\n").unwrap();

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
