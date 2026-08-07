use crate::error::{AppError, Result};
use crate::util::{revcomp, sanitize_filename_component};
use std::collections::HashSet;
use std::path::Path;

/// Demultiplex operating mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DemuxMode {
    FivePrime,
    FiveAndThreePrime,
    ThreePrimeOnly,
}

/// Compiled 5' barcode pattern.
#[derive(Debug, Clone)]
pub struct CompiledFivePrime {
    pub id: usize,
    pub raw: Vec<u8>,
    pub pattern_len: usize,
    pub informative_positions: Vec<usize>,
    pub expected_bases: Vec<u8>,
    pub umi_positions: Vec<usize>,
    pub sample_name: Option<String>,
    pub linked_three_prime: Vec<CompiledThreePrime>,
}

/// Compiled 3' barcode pattern (positions from read end).
#[derive(Debug, Clone)]
pub struct CompiledThreePrime {
    pub id: usize,
    pub raw: Vec<u8>,
    pub pattern_len: usize,
    /// Offsets from the 3' end: 0 = last base.
    pub informative_offsets_from_end: Vec<usize>,
    pub expected_bases: Vec<u8>,
    pub umi_offsets_from_end: Vec<usize>,
    pub sample_name: Option<String>,
}

/// Unique sample output key.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SampleKey {
    Named(String),
    FivePrime { barcode: String },
    Combined { five: String, three: String },
    ThreePrime { barcode: String },
    Unassigned,
}

impl SampleKey {
    pub fn label(&self) -> String {
        match self {
            SampleKey::Named(s) => sanitize_filename_component(s),
            SampleKey::FivePrime { barcode } => {
                format!("5bc_{}", sanitize_filename_component(barcode))
            }
            SampleKey::Combined { five, three } => format!(
                "5bc_{}_3bc_{}",
                sanitize_filename_component(five),
                sanitize_filename_component(three)
            ),
            SampleKey::ThreePrime { barcode } => {
                format!("3bc_{}", sanitize_filename_component(barcode))
            }
            SampleKey::Unassigned => "unassigned".to_string(),
        }
    }
}

/// Full barcode configuration after CSV parsing.
#[derive(Debug, Clone)]
pub struct BarcodeConfig {
    pub five_prime: Vec<CompiledFivePrime>,
    pub mode: DemuxMode,
    pub mismatches_5: usize,
    pub mismatches_3: usize,
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
    let upper = s.to_ascii_uppercase();
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

fn compile_five(raw: Vec<u8>, id: usize, sample_name: Option<String>) -> CompiledFivePrime {
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
    CompiledFivePrime {
        id,
        pattern_len: raw.len(),
        informative_positions,
        expected_bases,
        umi_positions,
        sample_name,
        linked_three_prime: Vec::new(),
        raw,
    }
}

fn compile_three(raw: Vec<u8>, id: usize, sample_name: Option<String>) -> CompiledThreePrime {
    let pattern_len = raw.len();
    let mut informative_offsets_from_end = Vec::new();
    let mut expected_bases = Vec::new();
    let mut umi_offsets_from_end = Vec::new();
    for (i, &b) in raw.iter().enumerate() {
        // offset from end: last base is 0
        let offset = pattern_len - 1 - i;
        if b == b'N' {
            umi_offsets_from_end.push(offset);
        } else {
            informative_offsets_from_end.push(offset);
            expected_bases.push(b);
        }
    }
    // Keep informative offsets sorted by position from 5'→3' of the pattern
    // (i.e. descending offset order as we walk from pattern start).
    // expected_bases already follows pattern 5'→3' order among non-N bases.
    CompiledThreePrime {
        id,
        raw,
        pattern_len,
        informative_offsets_from_end,
        expected_bases,
        umi_offsets_from_end,
        sample_name,
    }
}

/// Parse Ultraplex-compatible barcode CSV.
pub fn load_barcodes_csv<P: AsRef<Path>>(
    path: P,
    mismatches_5: usize,
    mismatches_3: usize,
    three_prime_only: bool,
) -> Result<BarcodeConfig> {
    let content = std::fs::read_to_string(path.as_ref())?;
    parse_barcodes_csv(&content, mismatches_5, mismatches_3, three_prime_only)
}

pub fn parse_barcodes_csv(
    content: &str,
    mismatches_5: usize,
    mismatches_3: usize,
    three_prime_only: bool,
) -> Result<BarcodeConfig> {
    let mut five_prime: Vec<CompiledFivePrime> = Vec::new();
    let mut sample_names_seen: HashSet<String> = HashSet::new();
    let mut five_info_len: Option<usize> = None;
    let mut has_any_three = false;
    let mut five_id = 0usize;
    let mut three_id = 0usize;

    for (line_no, raw_line) in content.lines().enumerate() {
        let line = raw_line.trim().replace(' ', "");
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split(',').collect();
        let first = parts[0];
        let (five_raw_str, five_sample) = split_bc_sample(first)?;
        let five_bytes = validate_bases(five_raw_str)?;

        let info_len = five_bytes.iter().filter(|&&b| b != b'N').count();
        if info_len == 0 {
            return Err(AppError::BarcodeConfig(format!(
                "line {}: barcode has no informative bases",
                line_no + 1
            )));
        }
        match five_info_len {
            None => five_info_len = Some(info_len),
            Some(n) if n != info_len => {
                return Err(AppError::BarcodeConfig(format!(
                    "line {}: 5' barcodes have inconsistent informative length ({n} vs {info_len})",
                    line_no + 1
                )));
            }
            _ => {}
        }

        // Collect linked 3' barcodes
        let mut linked = Vec::new();
        let mut three_info_len: Option<usize> = None;
        for col in parts.iter().skip(1) {
            if col.is_empty() {
                continue;
            }
            let (three_raw, three_sample) = split_bc_sample(col)?;
            let three_bytes = validate_bases(three_raw)?;
            let t_info = three_bytes.iter().filter(|&&b| b != b'N').count();
            if t_info == 0 {
                return Err(AppError::BarcodeConfig(format!(
                    "line {}: 3' barcode has no informative bases",
                    line_no + 1
                )));
            }
            match three_info_len {
                None => three_info_len = Some(t_info),
                Some(n) if n != t_info => {
                    return Err(AppError::BarcodeConfig(format!(
                        "line {}: linked 3' barcodes have inconsistent informative length",
                        line_no + 1
                    )));
                }
                _ => {}
            }
            if let Some(ref name) = three_sample {
                if !sample_names_seen.insert(name.clone()) {
                    return Err(AppError::BarcodeConfig(format!(
                        "duplicate sample name: {name}"
                    )));
                }
            }
            linked.push(compile_three(three_bytes, three_id, three_sample));
            three_id += 1;
        }

        if !linked.is_empty() {
            has_any_three = true;
            if five_sample.is_some() {
                return Err(AppError::BarcodeConfig(format!(
                    "line {}: cannot name 5' barcode when linked 3' barcodes are present",
                    line_no + 1
                )));
            }
        } else if let Some(ref name) = five_sample {
            if !sample_names_seen.insert(name.clone()) {
                return Err(AppError::BarcodeConfig(format!(
                    "duplicate sample name: {name}"
                )));
            }
        }

        let mut compiled = compile_five(five_bytes, five_id, five_sample);
        five_id += 1;
        compiled.linked_three_prime = linked;
        five_prime.push(compiled);
    }

    if five_prime.is_empty() {
        return Err(AppError::BarcodeConfig("no barcodes found in CSV".into()));
    }

    // Validate UMI position consistency for 5'
    check_five_n_positions(&five_prime)?;

    // For three_prime_only: reverse-complement patterns and require paired-end later
    let mode = if three_prime_only {
        if has_any_three {
            return Err(AppError::BarcodeConfig(
                "three_prime_only mode expects barcodes listed as 5' columns (no linked 3')".into(),
            ));
        }
        // RC each barcode
        for bc in &mut five_prime {
            let rc = revcomp(&bc.raw);
            let sample = bc.sample_name.clone();
            let id = bc.id;
            *bc = compile_five(rc, id, sample);
        }
        DemuxMode::ThreePrimeOnly
    } else if has_any_three {
        DemuxMode::FiveAndThreePrime
    } else {
        DemuxMode::FivePrime
    };

    let info5 = five_info_len.unwrap_or(0);
    if mismatches_5 > info5 {
        return Err(AppError::BarcodeConfig(format!(
            "mismatches-5 ({mismatches_5}) exceeds informative bases ({info5})"
        )));
    }
    if has_any_three {
        // check max 3' info across all
        let max3 = five_prime
            .iter()
            .flat_map(|f| f.linked_three_prime.iter())
            .map(|t| t.expected_bases.len())
            .max()
            .unwrap_or(0);
        if mismatches_3 > max3 {
            return Err(AppError::BarcodeConfig(format!(
                "mismatches-3 ({mismatches_3}) exceeds informative bases ({max3})"
            )));
        }
        for f in &five_prime {
            if !f.linked_three_prime.is_empty() {
                check_three_n_positions(&f.linked_three_prime)?;
            }
        }
    }

    Ok(BarcodeConfig {
        five_prime,
        mode,
        mismatches_5,
        mismatches_3,
    })
}

fn split_bc_sample(s: &str) -> Result<(&str, Option<String>)> {
    let colon_count = s.matches(':').count();
    if colon_count > 1 {
        return Err(AppError::BarcodeConfig(format!(
            "multiple colons in barcode field: {s}"
        )));
    }
    if let Some((bc, name)) = s.split_once(':') {
        let name = name.trim();
        if name.is_empty() {
            Ok((bc, None))
        } else {
            Ok((bc, Some(name.to_string())))
        }
    } else {
        Ok((s, None))
    }
}

fn check_five_n_positions(bcs: &[CompiledFivePrime]) -> Result<()> {
    if bcs.is_empty() {
        return Ok(());
    }
    // Ultraplex: first non-N position must be consistent
    let ref_pos = bcs[0]
        .informative_positions
        .first()
        .copied()
        .ok_or_else(|| AppError::BarcodeConfig("empty informative positions".into()))?;
    for bc in bcs {
        let pos = bc.informative_positions.first().copied().unwrap_or(0);
        if pos != ref_pos {
            return Err(AppError::BarcodeConfig(
                "UMI positions not consistent across 5' barcodes".into(),
            ));
        }
        // also require identical informative position sets for matching
        if bc.informative_positions != bcs[0].informative_positions {
            // Allow different UMI lengths only if informative positions match relative layout
            // Ultraplex requires same informative length; positions of non-N relative to first.
            // We require exact same informative positions for simplicity and correctness.
            return Err(AppError::BarcodeConfig(
                "5' barcode informative positions are not consistent".into(),
            ));
        }
    }
    Ok(())
}

fn check_three_n_positions(bcs: &[CompiledThreePrime]) -> Result<()> {
    if bcs.is_empty() {
        return Ok(());
    }
    // Ultraplex: distance from end of last non-N must be consistent
    let ref_off = bcs[0]
        .informative_offsets_from_end
        .iter()
        .min()
        .copied()
        .unwrap_or(0);
    for bc in bcs {
        let off = bc
            .informative_offsets_from_end
            .iter()
            .min()
            .copied()
            .unwrap_or(0);
        if off != ref_off {
            return Err(AppError::BarcodeConfig(
                "UMI positions not consistent across 3' barcodes".into(),
            ));
        }
    }
    Ok(())
}

/// Build sample key for a match.
pub fn sample_key_for(
    five: &CompiledFivePrime,
    three: Option<&CompiledThreePrime>,
    mode: DemuxMode,
) -> SampleKey {
    match mode {
        DemuxMode::ThreePrimeOnly => {
            if let Some(name) = &five.sample_name {
                SampleKey::Named(name.clone())
            } else {
                SampleKey::ThreePrime {
                    barcode: String::from_utf8_lossy(&five.raw).into_owned(),
                }
            }
        }
        _ => {
            if let Some(t) = three {
                if let Some(name) = &t.sample_name {
                    SampleKey::Named(name.clone())
                } else {
                    SampleKey::Combined {
                        five: String::from_utf8_lossy(&five.raw).into_owned(),
                        three: String::from_utf8_lossy(&t.raw).into_owned(),
                    }
                }
            } else if let Some(name) = &five.sample_name {
                SampleKey::Named(name.clone())
            } else {
                SampleKey::FivePrime {
                    barcode: String::from_utf8_lossy(&five.raw).into_owned(),
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_simple_5p() {
        let csv = "NNNATGNN:sampleA\nNNNCCGNN:sampleB\n";
        let cfg = parse_barcodes_csv(csv, 0, 0, false).unwrap();
        assert_eq!(cfg.five_prime.len(), 2);
        assert_eq!(cfg.mode, DemuxMode::FivePrime);
        assert_eq!(cfg.five_prime[0].sample_name.as_deref(), Some("sampleA"));
        assert_eq!(cfg.five_prime[0].expected_bases, b"ATG");
        assert_eq!(cfg.five_prime[0].umi_positions, vec![0, 1, 2, 6, 7]);
    }

    #[test]
    fn parse_linked_3p() {
        let csv = "NNNATGNN,ATG:s1,TCA:s2\nNNNCCGNN,\n";
        let cfg = parse_barcodes_csv(csv, 0, 0, false).unwrap();
        assert_eq!(cfg.mode, DemuxMode::FiveAndThreePrime);
        assert_eq!(cfg.five_prime[0].linked_three_prime.len(), 2);
        assert!(cfg.five_prime[1].linked_three_prime.is_empty());
    }

    #[test]
    fn reject_duplicate_sample() {
        let csv = "ATG:s1\nCCG:s1\n";
        assert!(parse_barcodes_csv(csv, 0, 0, false).is_err());
    }

    #[test]
    fn tso_parse() {
        let t = TsoPattern::parse("NNNNNIII").unwrap();
        assert_eq!(t.total_len, 8);
        assert_eq!(t.umi_positions.len(), 5);
        assert_eq!(t.ignored_positions.len(), 3);
    }
}
