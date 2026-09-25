use crate::error::{AppError, Result};
use crate::util::sanitize_filename_component;
use std::collections::{HashMap, HashSet};
use std::path::Path;

/// Demultiplex operating mode derived from the sample table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DemuxMode {
    /// Only Barcode1 present: match at R1 (or SE) 5' end.
    SingleBarcode,
    /// Barcode1 + Barcode2: PE → R1/R2 5'; SE → R1 5' + R1 3'.
    DualBarcode,
}

/// How dual-barcode PE reads may be oriented on R1/R2.
///
/// Amplicon libraries (I395 / I464) typically mix both orientations ~50/50
/// because inserts ligate in either direction relative to Illumina adapters.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OrientationMode {
    /// Try Barcode1@R1 + Barcode2@R2 and the swapped mates; pick the unique best.
    #[default]
    Both,
    /// Only Barcode1 at R1 5′ and Barcode2 at R2 5′.
    Canonical,
    /// Only Barcode2 at R1 5′ and Barcode1 at R2 5′.
    Swapped,
}

/// Compiled barcode pattern (supports optional `N` UMI bases).
#[derive(Debug, Clone)]
pub struct CompiledBarcode {
    pub raw: Vec<u8>,
    pub pattern_len: usize,
    pub informative_positions: Vec<usize>,
    pub expected_bases: Vec<u8>,
    pub umi_positions: Vec<usize>,
}

impl CompiledBarcode {
    pub fn compile(raw: Vec<u8>) -> Result<Self> {
        let mut informative_positions = Vec::new();
        let mut expected_bases = Vec::new();
        let mut umi_positions = Vec::new();
        for (i, &b) in raw.iter().enumerate() {
            if b == b'N' {
                umi_positions.push(i);
            } else {
                informative_positions.push(i);
                expected_bases.push(b);
            }
        }
        if expected_bases.is_empty() {
            return Err(AppError::BarcodeConfig(
                "barcode has no informative (non-N) bases".into(),
            ));
        }
        Ok(Self {
            pattern_len: raw.len(),
            informative_positions,
            expected_bases,
            umi_positions,
            raw,
        })
    }
}

/// One sample row from the barcode table.
#[derive(Debug, Clone)]
pub struct SampleEntry {
    pub id: usize,
    pub name: String,
    pub sanitized_label: String,
    pub barcode1: CompiledBarcode,
    pub barcode2: Option<CompiledBarcode>,
}

/// Compact 2-bit packing for 8-bp ACGT barcode:
/// A = 00, C = 01, G = 10, T = 11.
/// Returns None if length < 8 or contains non-ACGT bases (e.g. N).
#[inline(always)]
pub fn pack_8bp_2bit(seq: &[u8]) -> Option<u16> {
    if seq.len() < 8 {
        return None;
    }
    let mut code = 0u16;
    for &b in &seq[..8] {
        let v = match b {
            b'A' | b'a' => 0u16,
            b'C' | b'c' => 1u16,
            b'G' | b'g' => 2u16,
            b'T' | b't' => 3u16,
            _ => return None,
        };
        code = (code << 2) | v;
    }
    Some(code)
}

/// Exact fast path matcher for fixed 8-bp A/C/G/T dual barcodes with mismatch = 0.
#[derive(Debug, Clone)]
pub struct ExactDual8Matcher {
    /// Maps packed (BC1 << 16 | BC2) to sample index in samples list.
    pub table: HashMap<u32, usize>,
}

/// Unique sample output key.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SampleKey {
    Named(String),
    Unassigned,
}

impl SampleKey {
    pub fn label(&self) -> String {
        match self {
            SampleKey::Named(s) => sanitize_filename_component(s),
            SampleKey::Unassigned => "unassigned".to_string(),
        }
    }
}

/// Full barcode configuration after CSV parsing.
#[derive(Debug, Clone)]
pub struct BarcodeConfig {
    pub samples: Vec<SampleEntry>,
    pub mode: DemuxMode,
    /// Allowed mismatches for Barcode1.
    pub mismatches_1: usize,
    /// Allowed mismatches for Barcode2.
    pub mismatches_2: usize,
    /// Exact 8-bp dual-barcode fast path lookup table when applicable.
    pub fast_exact_8bp_pe: Option<ExactDual8Matcher>,
}

/// TSO pattern: N → UMI, I → trim only.
#[derive(Debug, Clone)]
pub struct TsoPattern {
    pub total_len: usize,
    pub umi_positions: Vec<usize>,
    pub ignored_positions: Vec<usize>,
}

impl TsoPattern {
    pub fn parse(s: &str) -> Result<Self> {
        let upper = s.to_ascii_uppercase();
        let mut umi_positions = Vec::new();
        let mut ignored_positions = Vec::new();
        for (i, c) in upper.bytes().enumerate() {
            match c {
                b'N' => umi_positions.push(i),
                b'I' => ignored_positions.push(i),
                other => {
                    return Err(AppError::BarcodeConfig(format!(
                        "invalid TSO base '{}', only N and I allowed",
                        other as char
                    )));
                }
            }
        }
        if upper.is_empty() {
            return Err(AppError::BarcodeConfig("empty TSO pattern".into()));
        }
        Ok(Self {
            total_len: upper.len(),
            umi_positions,
            ignored_positions,
        })
    }
}

fn validate_bases(s: &str) -> Result<Vec<u8>> {
    let upper = s.trim().to_ascii_uppercase();
    let bytes = upper.into_bytes();
    for &b in &bytes {
        if !matches!(b, b'A' | b'C' | b'G' | b'T' | b'N') {
            return Err(AppError::InvalidBase(format!(
                "invalid base '{}' in barcode (only A/C/G/T/N allowed)",
                b as char
            )));
        }
    }
    if bytes.is_empty() {
        return Err(AppError::BarcodeConfig("empty barcode".into()));
    }
    Ok(bytes)
}

fn normalize_header(h: &str) -> String {
    h.trim()
        .trim_start_matches('\u{feff}') // BOM
        .to_ascii_lowercase()
        .replace([' ', '-'], "_")
}

fn resolve_column(header: &str) -> Option<&'static str> {
    match normalize_header(header).as_str() {
        "samplenumber" | "sample_number" | "sample" | "sample_name" | "sampleid" | "sample_id"
        | "name" => Some("sample"),
        "barcode1" | "barcode_1" | "bc1" | "bc_1" | "index1" | "index_1" | "i7" => Some("barcode1"),
        "barcode2" | "barcode_2" | "bc2" | "bc_2" | "index2" | "index_2" | "i5" => Some("barcode2"),
        _ => None,
    }
}

/// Load SeqMux sample barcode table (CSV with header).
pub fn load_barcodes_csv<P: AsRef<Path>>(
    path: P,
    mismatches_1: usize,
    mismatches_2: usize,
) -> Result<BarcodeConfig> {
    let content = std::fs::read_to_string(path.as_ref()).map_err(|e| {
        AppError::BarcodeConfig(format!(
            "failed to read barcode table {}: {e}",
            path.as_ref().display()
        ))
    })?;
    parse_barcodes_csv(&content, mismatches_1, mismatches_2)
}

/// Parse SeqMux sample barcode CSV.
///
/// Required header columns (case-insensitive):
/// - `SampleNumber` (aliases: Sample, Sample_Name, Name, …)
/// - `Barcode1` (aliases: Barcode_1, BC1, Index1, i7, …)
///
/// Optional:
/// - `Barcode2` (aliases: Barcode_2, BC2, Index2, i5, …)
///
/// Extra columns (`PCR_product`, `rawdata1`, `F_primer`, …) are ignored.
pub fn parse_barcodes_csv(
    content: &str,
    mismatches_1: usize,
    mismatches_2: usize,
) -> Result<BarcodeConfig> {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .trim(csv::Trim::All)
        .from_reader(content.as_bytes());

    let headers = rdr.headers().map_err(|e| {
        AppError::BarcodeConfig(format!(
            "failed to read CSV header (SeqMux requires a header row with SampleNumber,Barcode1,Barcode2): {e}"
        ))
    })?;

    let mut col_sample: Option<usize> = None;
    let mut col_bc1: Option<usize> = None;
    let mut col_bc2: Option<usize> = None;

    for (i, h) in headers.iter().enumerate() {
        match resolve_column(h) {
            Some("sample") => col_sample = Some(i),
            Some("barcode1") => col_bc1 = Some(i),
            Some("barcode2") => col_bc2 = Some(i),
            _ => {}
        }
    }

    let col_sample = col_sample.ok_or_else(|| {
        AppError::BarcodeConfig(
            "missing required column SampleNumber (or Sample / Sample_Name / Name)".into(),
        )
    })?;
    let col_bc1 = col_bc1.ok_or_else(|| {
        AppError::BarcodeConfig(
            "missing required column Barcode1 (or Barcode_1 / BC1 / Index1)".into(),
        )
    })?;

    let mut samples = Vec::new();
    let mut names_seen: HashSet<String> = HashSet::new();
    let mut sanitized_seen: HashMap<String, String> = HashMap::new();
    let mut pair_seen: HashSet<(String, String)> = HashSet::new();
    let mut any_bc2 = false;
    let mut all_bc2 = true;

    for (row_idx, rec) in rdr.records().enumerate() {
        let rec = rec.map_err(|e| {
            AppError::BarcodeConfig(format!("CSV parse error at data row {}: {e}", row_idx + 1))
        })?;
        let line_no = row_idx + 2; // 1-based, accounting for header
        let id = row_idx;

        let name = rec.get(col_sample).unwrap_or("").trim().to_string();
        if name.is_empty() {
            return Err(AppError::BarcodeConfig(format!(
                "line {line_no}: empty SampleNumber"
            )));
        }
        if !names_seen.insert(name.clone()) {
            return Err(AppError::BarcodeConfig(format!(
                "line {line_no}: duplicate sample name '{name}'"
            )));
        }

        let label = sanitize_filename_component(&name);
        if label.eq_ignore_ascii_case("unassigned") {
            return Err(AppError::BarcodeConfig(format!(
                "line {line_no}: sample name '{name}' collides with reserved output name 'unassigned'"
            )));
        }
        if let Some(prev) = sanitized_seen.get(&label) {
            return Err(AppError::BarcodeConfig(format!(
                "line {line_no}: sample name '{name}' sanitizes to '{label}', which collides with earlier sample '{prev}'"
            )));
        }
        sanitized_seen.insert(label.clone(), name.clone());

        let bc1_raw = rec.get(col_bc1).unwrap_or("").trim();
        if bc1_raw.is_empty() {
            return Err(AppError::BarcodeConfig(format!(
                "line {line_no}: empty Barcode1 for sample '{name}'"
            )));
        }
        let barcode1 = CompiledBarcode::compile(validate_bases(bc1_raw)?)?;

        let barcode2 = if let Some(c2) = col_bc2 {
            let bc2_raw = rec.get(c2).unwrap_or("").trim();
            if bc2_raw.is_empty() {
                all_bc2 = false;
                None
            } else {
                any_bc2 = true;
                Some(CompiledBarcode::compile(validate_bases(bc2_raw)?)?)
            }
        } else {
            all_bc2 = false;
            None
        };

        let pair_key = (
            String::from_utf8_lossy(&barcode1.raw).into_owned(),
            barcode2
                .as_ref()
                .map(|b| String::from_utf8_lossy(&b.raw).into_owned())
                .unwrap_or_default(),
        );
        if !pair_seen.insert(pair_key) {
            return Err(AppError::BarcodeConfig(format!(
                "line {line_no}: duplicate barcode combination for sample '{name}'"
            )));
        }

        samples.push(SampleEntry {
            id,
            name,
            sanitized_label: label,
            barcode1,
            barcode2,
        });
    }

    if samples.is_empty() {
        return Err(AppError::BarcodeConfig(
            "no samples found in barcode table".into(),
        ));
    }

    if any_bc2 && !all_bc2 {
        return Err(AppError::BarcodeConfig(
            "Barcode2 must be present for all samples or none (mixed single/dual not supported)"
                .into(),
        ));
    }

    let mode = if all_bc2 && any_bc2 {
        DemuxMode::DualBarcode
    } else {
        DemuxMode::SingleBarcode
    };

    if mode == DemuxMode::DualBarcode {
        let set1: HashSet<&[u8]> = samples.iter().map(|s| s.barcode1.raw.as_slice()).collect();
        let set2: HashSet<&[u8]> = samples
            .iter()
            .filter_map(|s| s.barcode2.as_ref().map(|b| b.raw.as_slice()))
            .collect();
        let overlap: Vec<&[u8]> = set1.intersection(&set2).copied().collect();
        if !overlap.is_empty() {
            log::warn!(
                "{} barcode sequence(s) appear as both Barcode1 and Barcode2; \
                 --orientation both may produce extra ambiguous assignments",
                overlap.len()
            );
        }
    }

    let max_info1 = samples
        .iter()
        .map(|s| s.barcode1.expected_bases.len())
        .max()
        .unwrap_or(0);
    if mismatches_1 > max_info1 {
        return Err(AppError::BarcodeConfig(format!(
            "mismatches-1 ({mismatches_1}) exceeds longest Barcode1 informative length ({max_info1})"
        )));
    }
    if mode == DemuxMode::DualBarcode {
        let max_info2 = samples
            .iter()
            .filter_map(|s| s.barcode2.as_ref())
            .map(|b| b.expected_bases.len())
            .max()
            .unwrap_or(0);
        if mismatches_2 > max_info2 {
            return Err(AppError::BarcodeConfig(format!(
                "mismatches-2 ({mismatches_2}) exceeds longest Barcode2 informative length ({max_info2})"
            )));
        }
    }

    let fast_exact_8bp_pe =
        if mode == DemuxMode::DualBarcode && mismatches_1 == 0 && mismatches_2 == 0 {
            let mut table = HashMap::with_capacity(samples.len());
            let mut eligible = true;
            for (i, s) in samples.iter().enumerate() {
                let Some(bc2) = s.barcode2.as_ref() else {
                    eligible = false;
                    break;
                };
                if s.barcode1.pattern_len != 8 || bc2.pattern_len != 8 {
                    eligible = false;
                    break;
                }
                if !s.barcode1.umi_positions.is_empty() || !bc2.umi_positions.is_empty() {
                    eligible = false;
                    break;
                }
                let (Some(c1), Some(c2)) =
                    (pack_8bp_2bit(&s.barcode1.raw), pack_8bp_2bit(&bc2.raw))
                else {
                    eligible = false;
                    break;
                };
                let key = ((c1 as u32) << 16) | (c2 as u32);
                table.insert(key, i);
            }
            if eligible {
                Some(ExactDual8Matcher { table })
            } else {
                None
            }
        } else {
            None
        };

    Ok(BarcodeConfig {
        samples,
        mode,
        mismatches_1,
        mismatches_2,
        fast_exact_8bp_pe,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    const EXAMPLE: &str = "\
SampleNumber,Barcode1,Barcode2,PCR_product,rawdata1,rawdata2,library_round,F_primer,R_primer
I464469-A1,AAGTCCAA,GGAGTACT,atgc,E1.fq.gz,E2.fq.gz,A,fwd,rev
I464469-A2,GACCTGAA,GGACTTGG,agcg,E1.fq.gz,E2.fq.gz,A,fwd,rev
";

    #[test]
    fn parse_seqmux_table() {
        let cfg = parse_barcodes_csv(EXAMPLE, 0, 0).unwrap();
        assert_eq!(cfg.samples.len(), 2);
        assert_eq!(cfg.mode, DemuxMode::DualBarcode);
        assert_eq!(cfg.samples[0].name, "I464469-A1");
        assert_eq!(cfg.samples[0].barcode1.raw, b"AAGTCCAA");
        assert_eq!(cfg.samples[0].barcode2.as_ref().unwrap().raw, b"GGAGTACT");
        assert_eq!(cfg.samples[1].name, "I464469-A2");
    }

    #[test]
    fn parse_single_barcode() {
        let csv = "SampleNumber,Barcode1\ns1,ATGATGAT\ns2,CCGTAACG\n";
        let cfg = parse_barcodes_csv(csv, 0, 0).unwrap();
        assert_eq!(cfg.mode, DemuxMode::SingleBarcode);
        assert!(cfg.samples[0].barcode2.is_none());
    }

    #[test]
    fn reject_duplicate_sample() {
        let csv = "SampleNumber,Barcode1,Barcode2\ns1,AAAAAAAA,TTTTTTTT\ns1,CCCCCCCC,GGGGGGGG\n";
        assert!(parse_barcodes_csv(csv, 0, 0).is_err());
    }

    #[test]
    fn reject_missing_header() {
        let csv = "s1,AAAAAAAA,TTTTTTTT\n";
        assert!(parse_barcodes_csv(csv, 0, 0).is_err());
    }

    #[test]
    fn header_aliases() {
        let csv = "sample_name,bc1,bc2\nS1,ACGTACGT,TGCATGCA\n";
        let cfg = parse_barcodes_csv(csv, 0, 0).unwrap();
        assert_eq!(cfg.samples[0].name, "S1");
    }

    #[test]
    fn tso_parse() {
        let t = TsoPattern::parse("NNNNNIII").unwrap();
        assert_eq!(t.total_len, 8);
        assert_eq!(t.umi_positions.len(), 5);
        assert_eq!(t.ignored_positions.len(), 3);
    }

    #[test]
    fn parse_fixture_i464_table() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/I464-469erdai_barcode_and_name.csv");
        let cfg = load_barcodes_csv(&path, 0, 0).expect("parse fixture table");
        assert_eq!(cfg.samples.len(), 35);
        assert_eq!(cfg.mode, DemuxMode::DualBarcode);
        assert_eq!(cfg.samples[0].name, "I464469-A1");
        assert_eq!(cfg.samples[0].barcode1.raw, b"AAGTCCAA");
        assert_eq!(cfg.samples[0].barcode2.as_ref().unwrap().raw, b"GGAGTACT");
    }

    #[test]
    fn reject_sanitized_sample_name_collision() {
        // "sample:1" and "sample/1" both sanitize to "sample_1"
        let csv = "SampleNumber,Barcode1,Barcode2\nsample:1,AAGTCCAA,GGAGTACT\nsample/1,GACCTGAA,GGACTTGG\n";
        let err = parse_barcodes_csv(csv, 0, 0).unwrap_err();
        assert!(
            err.to_string().contains("collides"),
            "expected collision error, got: {err}"
        );
    }

    #[test]
    fn reject_unassigned_sample_name() {
        let csv = "SampleNumber,Barcode1,Barcode2\nUnassigned,AAGTCCAA,GGAGTACT\n";
        let err = parse_barcodes_csv(csv, 0, 0).unwrap_err();
        assert!(
            err.to_string().contains("unassigned"),
            "expected unassigned error, got: {err}"
        );
    }
}
