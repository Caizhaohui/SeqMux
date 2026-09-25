use crate::barcode::{
    assign_paired, assign_single, extract_umi_3p, extract_umi_5p, BarcodeConfig, DemuxMode,
    MatchResult, OrientationMode, TsoPattern,
};
use crate::fastq::{OwnedFastqRecord, ReadOrPair};
use crate::output::{Mate, OutputKey, UNASSIGNED_SAMPLE_ID};
use crate::stats::ChunkStats;
use crate::trim::{nextseq_trim_index, quality_trim_bounds, CompiledAdapter};
use std::collections::HashMap;

/// Processing parameters shared by workers.
#[derive(Debug, Clone)]
pub struct ProcessConfig {
    pub barcodes: BarcodeConfig,
    pub keep_barcodes: bool,
    pub discard_unassigned: bool,
    pub min_length: usize,
    pub quality_cutoff_5: u8,
    pub quality_cutoff_3: u8,
    pub nextseq: bool,
    pub adapter_r1: Option<CompiledAdapter>,
    pub adapter_r2: Option<CompiledAdapter>,
    pub adapter_error_rate: f32,
    pub min_adapter_overlap: usize,
    pub tso: Option<TsoPattern>,
    pub phred_offset: u8,
    /// Dual-barcode PE orientation policy.
    pub orientation: OrientationMode,
    /// When a swapped orientation matches, swap mates so output R1 has Barcode1.
    pub canonicalize: bool,
    /// Count assignments but do not write FASTQ (real-data QC / smoke tests).
    pub counts_only: bool,
}

#[derive(Debug)]
pub struct InputChunk {
    pub id: u64,
    pub records: Vec<ReadOrPair>,
}

#[derive(Debug)]
pub struct ProcessedChunk {
    pub id: u64,
    pub outputs: HashMap<OutputKey, Vec<u8>>,
    pub stats: ChunkStats,
}

pub fn process_chunk(chunk: InputChunk, cfg: &ProcessConfig) -> ProcessedChunk {
    let mut outputs: HashMap<OutputKey, Vec<u8>> = HashMap::new();
    let mut stats = ChunkStats::new(cfg.barcodes.samples.len());

    for item in chunk.records {
        match item {
            ReadOrPair::Single(rec) => {
                process_single(rec, cfg, &mut outputs, &mut stats);
            }
            ReadOrPair::Pair(pair) => {
                process_pair(pair.r1, pair.r2, cfg, &mut outputs, &mut stats);
            }
        }
    }

    ProcessedChunk {
        id: chunk.id,
        outputs,
        stats,
    }
}

fn quality_trim_record(rec: &mut OwnedFastqRecord, cfg: &ProcessConfig) -> bool {
    let orig = rec.len();
    if cfg.nextseq {
        let stop = nextseq_trim_index(
            &rec.sequence,
            &rec.qualities,
            cfg.quality_cutoff_3,
            cfg.phred_offset,
        );
        let trimmed = stop < orig;
        rec.trim_to(0, stop);
        trimmed
    } else {
        let (start, end) = quality_trim_bounds(
            &rec.qualities,
            cfg.quality_cutoff_5,
            cfg.quality_cutoff_3,
            cfg.phred_offset,
        );
        let trimmed = start > 0 || end < orig;
        rec.trim_to(start, end);
        trimmed
    }
}

fn adapter_trim_record(rec: &mut OwnedFastqRecord, adapter: Option<&CompiledAdapter>) -> bool {
    let Some(adapter) = adapter else {
        return false;
    };
    let (new_end, trimmed, _) = adapter.trim_3p(&rec.sequence);
    if trimmed {
        rec.trim_to(0, new_end);
    }
    trimmed
}

fn process_single(
    mut rec: OwnedFastqRecord,
    cfg: &ProcessConfig,
    outputs: &mut HashMap<OutputKey, Vec<u8>>,
    stats: &mut ChunkStats,
) {
    stats.total_reads += 1;
    // Match on the original 5′ sequence before quality/adapter trim.
    let result = assign_single(&rec.sequence, &cfg.barcodes);
    let (sample_idx, umi, ambiguous) = apply_assignment_single(&mut rec, result, cfg);

    if ambiguous {
        stats.ambiguous += 1;
    }
    if let Some(ref tso) = cfg.tso {
        apply_tso(&mut rec, tso);
    }
    if quality_trim_record(&mut rec, cfg) {
        stats.quality_trimmed += 1;
    }
    if adapter_trim_record(&mut rec, cfg.adapter_r1.as_ref()) {
        stats.adapter_trimmed += 1;
    }
    if rec.len() < cfg.min_length {
        stats.too_short += 1;
        return;
    }
    if !umi.is_empty() {
        rec.append_umi(&umi);
    }

    if sample_idx.is_none() && cfg.discard_unassigned {
        stats.record_sample_opt(sample_idx);
        return;
    }
    stats.record_sample_opt(sample_idx);
    if cfg.counts_only {
        return;
    }
    let key = OutputKey {
        sample_id: sample_idx.unwrap_or(UNASSIGNED_SAMPLE_ID),
        mate: Mate::Single,
    };
    rec.append_fastq_to(outputs.entry(key).or_default());
}

fn process_pair(
    mut r1: OwnedFastqRecord,
    mut r2: OwnedFastqRecord,
    cfg: &ProcessConfig,
    outputs: &mut HashMap<OutputKey, Vec<u8>>,
    stats: &mut ChunkStats,
) {
    stats.total_reads += 1;
    // Match on original 5′ bases before quality/adapter trim.
    let result = assign_paired(&r1.sequence, &r2.sequence, &cfg.barcodes, cfg.orientation);
    match result {
        MatchResult::Match { swapped: true, .. } => stats.orientation_swapped += 1,
        MatchResult::Match { swapped: false, .. } => stats.orientation_canonical += 1,
        _ => {}
    }
    let (sample_idx, umi, ambiguous) = apply_assignment_paired(&mut r1, &mut r2, result, cfg);

    if ambiguous {
        stats.ambiguous += 1;
    }
    if let Some(ref tso) = cfg.tso {
        apply_tso(&mut r1, tso);
    }
    let q1 = quality_trim_record(&mut r1, cfg);
    let q2 = quality_trim_record(&mut r2, cfg);
    if q1 || q2 {
        stats.quality_trimmed += 1;
    }
    let a1 = adapter_trim_record(&mut r1, cfg.adapter_r1.as_ref());
    let a2 = adapter_trim_record(&mut r2, cfg.adapter_r2.as_ref());
    if a1 || a2 {
        stats.adapter_trimmed += 1;
    }
    if r1.len() < cfg.min_length || r2.len() < cfg.min_length {
        stats.too_short += 1;
        return;
    }
    if !umi.is_empty() {
        r1.append_umi(&umi);
        r2.append_umi(&umi);
    }

    if sample_idx.is_none() && cfg.discard_unassigned {
        stats.record_sample_opt(sample_idx);
        return;
    }
    stats.record_sample_opt(sample_idx);
    if cfg.counts_only {
        return;
    }
    let sid = sample_idx.unwrap_or(UNASSIGNED_SAMPLE_ID);
    let k1 = OutputKey {
        sample_id: sid,
        mate: Mate::R1,
    };
    let k2 = OutputKey {
        sample_id: sid,
        mate: Mate::R2,
    };
    r1.append_fastq_to(outputs.entry(k1).or_default());
    r2.append_fastq_to(outputs.entry(k2).or_default());
}

fn apply_assignment_single(
    rec: &mut OwnedFastqRecord,
    result: MatchResult,
    cfg: &ProcessConfig,
) -> (Option<usize>, Vec<u8>, bool) {
    match result {
        MatchResult::NoMatch => (None, Vec::new(), false),
        MatchResult::Ambiguous { .. } => (None, Vec::new(), true),
        MatchResult::Match { sample_index, .. } => {
            let sample = &cfg.barcodes.samples[sample_index];
            let mut umi = extract_umi_5p(&rec.sequence, &sample.barcode1);
            let bc1_len = sample.barcode1.pattern_len;
            let bc2_len = sample.barcode2.as_ref().map(|b| b.pattern_len).unwrap_or(0);

            if matches!(cfg.barcodes.mode, DemuxMode::DualBarcode) {
                if let Some(bc2) = sample.barcode2.as_ref() {
                    umi.extend_from_slice(&extract_umi_3p(&rec.sequence, bc2));
                }
            }

            if !cfg.keep_barcodes {
                // Trim 3' barcode first so 5' indices stay valid
                if matches!(cfg.barcodes.mode, DemuxMode::DualBarcode)
                    && bc2_len > 0
                    && rec.len() >= bc1_len + bc2_len
                {
                    rec.trim_back(bc2_len);
                }
                if rec.len() >= bc1_len {
                    rec.trim_front(bc1_len);
                }
            }
            (Some(sample_index), umi, false)
        }
    }
}

fn apply_assignment_paired(
    r1: &mut OwnedFastqRecord,
    r2: &mut OwnedFastqRecord,
    result: MatchResult,
    cfg: &ProcessConfig,
) -> (Option<usize>, Vec<u8>, bool) {
    match result {
        MatchResult::NoMatch => (None, Vec::new(), false),
        MatchResult::Ambiguous { .. } => (None, Vec::new(), true),
        MatchResult::Match {
            sample_index,
            swapped,
            ..
        } => {
            let sample = &cfg.barcodes.samples[sample_index];
            let (umi, trim_swapped) = if swapped {
                let mut umi = extract_umi_5p(&r2.sequence, &sample.barcode1);
                if let Some(bc2) = sample.barcode2.as_ref() {
                    umi.extend_from_slice(&extract_umi_5p(&r1.sequence, bc2));
                }
                if cfg.canonicalize {
                    std::mem::swap(r1, r2);
                    (umi, false)
                } else {
                    (umi, true)
                }
            } else {
                let mut umi = extract_umi_5p(&r1.sequence, &sample.barcode1);
                if let Some(bc2) = sample.barcode2.as_ref() {
                    umi.extend_from_slice(&extract_umi_5p(&r2.sequence, bc2));
                }
                (umi, false)
            };

            if !cfg.keep_barcodes {
                if trim_swapped {
                    if let Some(bc2) = sample.barcode2.as_ref() {
                        let l2 = bc2.pattern_len;
                        if r1.len() >= l2 {
                            r1.trim_front(l2);
                        }
                    }
                    let l1 = sample.barcode1.pattern_len;
                    if r2.len() >= l1 {
                        r2.trim_front(l1);
                    }
                } else {
                    let l1 = sample.barcode1.pattern_len;
                    if r1.len() >= l1 {
                        r1.trim_front(l1);
                    }
                    if let Some(bc2) = sample.barcode2.as_ref() {
                        let l2 = bc2.pattern_len;
                        if r2.len() >= l2 {
                            r2.trim_front(l2);
                        }
                    }
                }
            }
            (Some(sample_index), umi, false)
        }
    }
}

fn apply_tso(rec: &mut OwnedFastqRecord, tso: &TsoPattern) {
    if rec.len() < tso.total_len {
        return;
    }
    let umi: Vec<u8> = tso
        .umi_positions
        .iter()
        .filter_map(|&p| rec.sequence.get(p).copied())
        .map(|b| b.to_ascii_uppercase())
        .collect();
    rec.trim_front(tso.total_len);
    if !umi.is_empty() {
        rec.append_umi(&umi);
    }
}
