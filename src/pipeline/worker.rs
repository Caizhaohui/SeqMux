use crate::barcode::{
    extract_five_umi, extract_three_umi, match_five_prime, match_three_prime, sample_key_for,
    BarcodeConfig, DemuxMode, MatchResult, SampleKey, TsoPattern,
};
use crate::fastq::{OwnedFastqRecord, ReadOrPair};
use crate::output::{Mate, OutputKey};
use crate::stats::ChunkStats;
use crate::trim::{nextseq_trim_index, quality_trim_bounds, trim_3p_adapter};
use crate::util::revcomp;
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
    pub adapter_r1: Option<Vec<u8>>,
    pub adapter_r2: Option<Vec<u8>>,
    pub adapter_error_rate: f32,
    pub min_adapter_overlap: usize,
    pub min_adapter_trim_for_3p: usize,
    pub tso: Option<TsoPattern>,
    pub phred_offset: u8,
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
    let mut stats = ChunkStats::default();

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

fn process_single(
    mut rec: OwnedFastqRecord,
    cfg: &ProcessConfig,
    outputs: &mut HashMap<OutputKey, Vec<u8>>,
    stats: &mut ChunkStats,
) {
    stats.total_reads += 1;

    // 1. Quality trim
    let orig_len = rec.len();
    if cfg.nextseq {
        let stop = nextseq_trim_index(
            &rec.sequence,
            &rec.qualities,
            cfg.quality_cutoff_3,
            cfg.phred_offset,
        );
        if stop < orig_len {
            stats.quality_trimmed += 1;
        }
        rec.trim_to(0, stop);
    } else {
        let (start, end) = quality_trim_bounds(
            &rec.qualities,
            cfg.quality_cutoff_5,
            cfg.quality_cutoff_3,
            cfg.phred_offset,
        );
        if start > 0 || end < orig_len {
            stats.quality_trimmed += 1;
        }
        rec.trim_to(start, end);
    }

    // 2. Adapter trim
    let mut adapter_overlap = 0usize;
    if let Some(ref adapter) = cfg.adapter_r1 {
        let (new_end, trimmed, ov) = trim_3p_adapter(
            &rec.sequence,
            adapter,
            cfg.adapter_error_rate,
            cfg.min_adapter_overlap,
        );
        if trimmed {
            stats.adapter_trimmed += 1;
            adapter_overlap = ov;
            rec.trim_to(0, new_end);
        }
    }

    // 3. Barcode demux
    let (sample, umi, ambiguous) = demux_read(
        &mut rec,
        cfg,
        adapter_overlap,
        true, // is_r1 / single
    );

    if ambiguous {
        stats.ambiguous += 1;
    }

    // 4. TSO (optional, after 5' demux for PE primarily; apply if set)
    if let Some(ref tso) = cfg.tso {
        apply_tso(&mut rec, tso);
    }

    // 5. Length filter
    if rec.len() < cfg.min_length {
        stats.too_short += 1;
        return;
    }

    // 6. UMI header
    if !umi.is_empty() {
        rec.append_umi(&umi);
    } else {
        rec.ensure_rbc_tag();
    }

    if matches!(sample, SampleKey::Unassigned) && cfg.discard_unassigned {
        stats.record_sample(&sample);
        return;
    }

    stats.record_sample(&sample);
    let key = OutputKey {
        sample,
        mate: Mate::Single,
    };
    outputs
        .entry(key)
        .or_default()
        .extend_from_slice(&rec.to_fastq_bytes());
}

fn process_pair(
    mut r1: OwnedFastqRecord,
    mut r2: OwnedFastqRecord,
    cfg: &ProcessConfig,
    outputs: &mut HashMap<OutputKey, Vec<u8>>,
    stats: &mut ChunkStats,
) {
    stats.total_reads += 1;

    // Quality trim both
    for rec in [&mut r1, &mut r2] {
        let orig = rec.len();
        if cfg.nextseq {
            let stop = nextseq_trim_index(
                &rec.sequence,
                &rec.qualities,
                cfg.quality_cutoff_3,
                cfg.phred_offset,
            );
            if stop < orig {
                stats.quality_trimmed += 1;
            }
            rec.trim_to(0, stop);
        } else {
            let (start, end) = quality_trim_bounds(
                &rec.qualities,
                cfg.quality_cutoff_5,
                cfg.quality_cutoff_3,
                cfg.phred_offset,
            );
            if start > 0 || end < orig {
                stats.quality_trimmed += 1;
            }
            rec.trim_to(start, end);
        }
    }

    // Adapter trim
    if let Some(ref adapter) = cfg.adapter_r1 {
        let (new_end, trimmed, _) = trim_3p_adapter(
            &r1.sequence,
            adapter,
            cfg.adapter_error_rate,
            cfg.min_adapter_overlap,
        );
        if trimmed {
            stats.adapter_trimmed += 1;
            r1.trim_to(0, new_end);
        }
    }
    if let Some(ref adapter) = cfg.adapter_r2 {
        let (new_end, trimmed, _) = trim_3p_adapter(
            &r2.sequence,
            adapter,
            cfg.adapter_error_rate,
            cfg.min_adapter_overlap,
        );
        if trimmed {
            stats.adapter_trimmed += 1;
            r2.trim_to(0, new_end);
        }
    }

    let (sample, umi, ambiguous) = match cfg.barcodes.mode {
        DemuxMode::ThreePrimeOnly => {
            // Barcodes are already RC'd; match on R2 (which has the RC barcode at 5')
            demux_three_prime_only(&mut r1, &mut r2, cfg)
        }
        _ => demux_paired(&mut r1, &mut r2, cfg),
    };

    if ambiguous {
        stats.ambiguous += 1;
    }

    if let Some(ref tso) = cfg.tso {
        apply_tso(&mut r1, tso);
    }

    if r1.len() < cfg.min_length || r2.len() < cfg.min_length {
        stats.too_short += 1;
        return;
    }

    if !umi.is_empty() {
        r1.append_umi(&umi);
        r2.append_umi(&umi);
    } else {
        r1.ensure_rbc_tag();
        r2.ensure_rbc_tag();
    }

    if matches!(sample, SampleKey::Unassigned) && cfg.discard_unassigned {
        stats.record_sample(&sample);
        return;
    }

    stats.record_sample(&sample);
    let k1 = OutputKey {
        sample: sample.clone(),
        mate: Mate::R1,
    };
    let k2 = OutputKey {
        sample,
        mate: Mate::R2,
    };
    outputs
        .entry(k1)
        .or_default()
        .extend_from_slice(&r1.to_fastq_bytes());
    outputs
        .entry(k2)
        .or_default()
        .extend_from_slice(&r2.to_fastq_bytes());
}

fn demux_read(
    rec: &mut OwnedFastqRecord,
    cfg: &ProcessConfig,
    adapter_overlap: usize,
    _is_single: bool,
) -> (SampleKey, Vec<u8>, bool) {
    let bcs = &cfg.barcodes;
    let five_match = match_five_prime(&rec.sequence, &bcs.five_prime, bcs.mismatches_5);

    match five_match {
        MatchResult::NoMatch => (SampleKey::Unassigned, Vec::new(), false),
        MatchResult::Ambiguous { .. } => (SampleKey::Unassigned, Vec::new(), true),
        MatchResult::Match { index, .. } => {
            let five = &bcs.five_prime[index];
            let mut umi = extract_five_umi(&rec.sequence, five);

            if five.pattern_len > rec.len() {
                return (SampleKey::Unassigned, Vec::new(), false);
            }

            if !cfg.keep_barcodes {
                rec.trim_front(five.pattern_len);
            }

            if five.linked_three_prime.is_empty() {
                let key = sample_key_for(five, None, bcs.mode);
                return (key, umi, false);
            }

            // Linked 3' barcode — for single-end, require min adapter trim
            if adapter_overlap < cfg.min_adapter_trim_for_3p {
                // 5' matched but cannot safely assign 3'
                return (SampleKey::Unassigned, umi, false);
            }

            let three_match =
                match_three_prime(&rec.sequence, &five.linked_three_prime, bcs.mismatches_3);
            match three_match {
                MatchResult::Match { index: ti, .. } => {
                    let three = &five.linked_three_prime[ti];
                    let umi3 = extract_three_umi(&rec.sequence, three);
                    umi.extend_from_slice(&umi3);
                    if !cfg.keep_barcodes && rec.len() >= three.pattern_len {
                        rec.trim_back(three.pattern_len);
                    }
                    let key = sample_key_for(five, Some(three), bcs.mode);
                    (key, umi, false)
                }
                MatchResult::Ambiguous { .. } => (SampleKey::Unassigned, umi, true),
                MatchResult::NoMatch => (SampleKey::Unassigned, umi, false),
            }
        }
    }
}

fn demux_paired(
    r1: &mut OwnedFastqRecord,
    r2: &mut OwnedFastqRecord,
    cfg: &ProcessConfig,
) -> (SampleKey, Vec<u8>, bool) {
    let bcs = &cfg.barcodes;
    let five_match = match_five_prime(&r1.sequence, &bcs.five_prime, bcs.mismatches_5);

    match five_match {
        MatchResult::NoMatch => (SampleKey::Unassigned, Vec::new(), false),
        MatchResult::Ambiguous { .. } => (SampleKey::Unassigned, Vec::new(), true),
        MatchResult::Match { index, .. } => {
            let five = &bcs.five_prime[index];
            if five.pattern_len > r1.len() {
                return (SampleKey::Unassigned, Vec::new(), false);
            }
            let mut umi = extract_five_umi(&r1.sequence, five);

            // Mate trim: remove RC of barcode (+ optional adapter context) from R2
            // Ultraplex removes the reverse complement of the 5' barcode from R2 3' end
            if !cfg.keep_barcodes {
                let rc = revcomp(&five.raw);
                // Try exact RC suffix on R2
                if r2.sequence.len() >= rc.len() {
                    let start = r2.sequence.len() - rc.len();
                    let suffix: Vec<u8> = r2.sequence[start..]
                        .iter()
                        .map(|b| b.to_ascii_uppercase())
                        .collect();
                    if suffix == rc {
                        r2.trim_back(rc.len());
                    }
                }
                r1.trim_front(five.pattern_len);
            }

            if five.linked_three_prime.is_empty() {
                let key = sample_key_for(five, None, bcs.mode);
                return (key, umi, false);
            }

            // 3' barcode on R1
            let three_match =
                match_three_prime(&r1.sequence, &five.linked_three_prime, bcs.mismatches_3);
            match three_match {
                MatchResult::Match { index: ti, .. } => {
                    let three = &five.linked_three_prime[ti];
                    let umi3 = extract_three_umi(&r1.sequence, three);
                    umi.extend_from_slice(&umi3);
                    if !cfg.keep_barcodes && r1.len() >= three.pattern_len {
                        r1.trim_back(three.pattern_len);
                    }
                    // Also try RC of 3' barcode at R2 5' end
                    if !cfg.keep_barcodes {
                        let rc3 = revcomp(&three.raw);
                        if r2.sequence.len() >= rc3.len() {
                            let prefix: Vec<u8> = r2.sequence[..rc3.len()]
                                .iter()
                                .map(|b| b.to_ascii_uppercase())
                                .collect();
                            if prefix == rc3 {
                                r2.trim_front(rc3.len());
                            }
                        }
                    }
                    let key = sample_key_for(five, Some(three), bcs.mode);
                    (key, umi, false)
                }
                MatchResult::Ambiguous { .. } => (SampleKey::Unassigned, umi, true),
                MatchResult::NoMatch => (SampleKey::Unassigned, umi, false),
            }
        }
    }
}

fn demux_three_prime_only(
    r1: &mut OwnedFastqRecord,
    r2: &mut OwnedFastqRecord,
    cfg: &ProcessConfig,
) -> (SampleKey, Vec<u8>, bool) {
    // Patterns already reverse-complemented at load time; match against R2 5' end
    let bcs = &cfg.barcodes;
    let five_match = match_five_prime(&r2.sequence, &bcs.five_prime, bcs.mismatches_5);
    match five_match {
        MatchResult::NoMatch => (SampleKey::Unassigned, Vec::new(), false),
        MatchResult::Ambiguous { .. } => (SampleKey::Unassigned, Vec::new(), true),
        MatchResult::Match { index, .. } => {
            let five = &bcs.five_prime[index];
            if five.pattern_len > r2.len() {
                return (SampleKey::Unassigned, Vec::new(), false);
            }
            let umi = extract_five_umi(&r2.sequence, five);
            if !cfg.keep_barcodes {
                r2.trim_front(five.pattern_len);
                // Original (non-RC) pattern was stored... we RC'd at load, so original is RC of current
                let orig = revcomp(&five.raw);
                if r1.sequence.len() >= orig.len() {
                    let start = r1.sequence.len() - orig.len();
                    let suffix: Vec<u8> = r1.sequence[start..]
                        .iter()
                        .map(|b| b.to_ascii_uppercase())
                        .collect();
                    if suffix == orig {
                        r1.trim_back(orig.len());
                    }
                }
            }
            let key = sample_key_for(five, None, DemuxMode::ThreePrimeOnly);
            (key, umi, false)
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
