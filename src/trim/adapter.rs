/// 3' adapter match result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdapterMatch {
    pub read_start: usize,
    pub read_end: usize,
    pub adapter_start: usize,
    pub adapter_end: usize,
    pub errors: usize,
    pub overlap: usize,
}

/// Open-addressing table for packed 8-mer → adapter-offset lookup (stride-1 complete).
const SEED_TABLE_SIZE: usize = 64;
const SEED_TABLE_MASK: usize = SEED_TABLE_SIZE - 1;
const EMPTY_OFF: u8 = 0xFF;
const SEED_OVERFLOW: usize = 8;
/// Fixed short 4-mer seeds covering adapter offsets that can participate in overlap < 16.
const SHORT4_CAP: usize = 16;

#[derive(Debug, Clone)]
struct SeedLookup {
    keys: [u64; SEED_TABLE_SIZE],
    offs: [u8; SEED_TABLE_SIZE],
    overflow: [(u64, u8); SEED_OVERFLOW],
    overflow_len: u8,
    n_seeds: u16,
}

impl SeedLookup {
    fn empty() -> Self {
        Self {
            keys: [0; SEED_TABLE_SIZE],
            offs: [EMPTY_OFF; SEED_TABLE_SIZE],
            overflow: [(0, 0); SEED_OVERFLOW],
            overflow_len: 0,
            n_seeds: 0,
        }
    }

    #[inline]
    fn hash(key: u64) -> usize {
        key.wrapping_mul(0x9E37_79B9_7F4A_7C15) as usize & SEED_TABLE_MASK
    }

    fn insert(&mut self, key: u64, offset: usize) {
        let off = offset as u8;
        let mut i = Self::hash(key);
        for _ in 0..SEED_TABLE_SIZE {
            if self.offs[i] == EMPTY_OFF {
                self.keys[i] = key;
                self.offs[i] = off;
                self.n_seeds += 1;
                return;
            }
            if self.keys[i] == key {
                if (self.overflow_len as usize) < SEED_OVERFLOW {
                    self.overflow[self.overflow_len as usize] = (key, off);
                    self.overflow_len += 1;
                    self.n_seeds += 1;
                }
                return;
            }
            i = (i + 1) & SEED_TABLE_MASK;
        }
        if (self.overflow_len as usize) < SEED_OVERFLOW {
            self.overflow[self.overflow_len as usize] = (key, off);
            self.overflow_len += 1;
            self.n_seeds += 1;
        }
    }

    #[inline]
    fn lookup(&self, key: u64, out: &mut [u8]) -> usize {
        if self.n_seeds == 0 || out.is_empty() {
            return 0;
        }
        let mut n = 0usize;
        let mut i = Self::hash(key);
        for _ in 0..SEED_TABLE_SIZE {
            if self.offs[i] == EMPTY_OFF {
                break;
            }
            if self.keys[i] == key {
                out[n] = self.offs[i];
                n += 1;
                if n == out.len() {
                    return n;
                }
            }
            i = (i + 1) & SEED_TABLE_MASK;
        }
        for j in 0..self.overflow_len as usize {
            if self.overflow[j].0 == key {
                out[n] = self.overflow[j].1;
                n += 1;
                if n == out.len() {
                    return n;
                }
            }
        }
        n
    }

    #[inline]
    fn is_empty(&self) -> bool {
        self.n_seeds == 0
    }
}

/// Precompiled adapter representation for fast zero-allocation 3' trimming.
#[derive(Debug, Clone)]
pub struct CompiledAdapter {
    pub seq: Vec<u8>,
    /// Stride-1 packed 8-mer → offset lookup (all discoverable seeds retained).
    seed_lookup: SeedLookup,
    /// 4-mer seeds for short-overlap indel rescue (Path 3), offsets 0..15.
    short4: [(u32, u8); SHORT4_CAP],
    short4_len: u8,
    pub max_error_rate: f32,
    pub min_overlap: usize,
}

impl CompiledAdapter {
    /// Precompile an adapter sequence and indexing structures.
    pub fn new(adapter: &[u8], max_error_rate: f32, min_overlap: usize) -> Self {
        let seq: Vec<u8> = adapter.iter().map(|b| b.to_ascii_uppercase()).collect();
        let mut seed_lookup = SeedLookup::empty();
        if seq.len() >= 8 {
            // Stride-1: every adapter 8-mer offset required by the correctness proof.
            for offset in 0..=seq.len() - 8 {
                let s = u64::from_ne_bytes(seq[offset..offset + 8].try_into().unwrap());
                seed_lookup.insert(s, offset);
            }
        }
        let mut short4 = [(0u32, 0u8); SHORT4_CAP];
        let mut short4_len = 0u8;
        if seq.len() >= 4 {
            let max_off = (seq.len() - 4).min(SHORT4_CAP - 1);
            for offset in 0..=max_off {
                let s = u32::from_ne_bytes(seq[offset..offset + 4].try_into().unwrap());
                short4[short4_len as usize] = (s, offset as u8);
                short4_len += 1;
            }
        }
        Self {
            seq,
            seed_lookup,
            short4,
            short4_len,
            max_error_rate,
            min_overlap,
        }
    }

    /// Access the uppercase adapter sequence bytes.
    #[inline]
    pub fn sequence(&self) -> &[u8] {
        &self.seq
    }

    /// Trim 3' adapter from read if found; returns (new_seq_end, was_trimmed, overlap).
    #[inline]
    pub fn trim_3p(&self, read: &[u8]) -> (usize, bool, usize) {
        if self.seq.is_empty() {
            return (read.len(), false, 0);
        }
        match self.find_3p(read) {
            Some(m) => (m.read_start, true, m.overlap),
            None => (read.len(), false, 0),
        }
    }

    /// Find best 3' adapter match in the read using fast seed-and-extend with stack DP.
    #[inline]
    pub fn find_3p(&self, read: &[u8]) -> Option<AdapterMatch> {
        find_3p_adapter_fast(read, self)
    }
}

/// Compare two candidate matches according to priority:
/// 1. More overlap
/// 2. Fewer errors
/// 3. Earlier read_start
#[inline]
pub fn prefer(a: AdapterMatch, b: AdapterMatch) -> AdapterMatch {
    if b.overlap > a.overlap {
        b
    } else if b.overlap < a.overlap {
        a
    } else if b.errors < a.errors {
        b
    } else if b.errors > a.errors {
        a
    } else if b.read_start < a.read_start {
        b
    } else {
        a
    }
}

/// Generic dynamic programming 3' adapter matcher (reference/fallback implementation).
pub fn find_3p_adapter_generic(
    read: &[u8],
    adapter: &[u8],
    max_error_rate: f32,
    min_overlap: usize,
) -> Option<AdapterMatch> {
    if adapter.is_empty() || read.is_empty() || min_overlap == 0 {
        return None;
    }
    let adapter: Vec<u8> = adapter.iter().map(|b| b.to_ascii_uppercase()).collect();
    let read_u: Vec<u8> = read.iter().map(|b| b.to_ascii_uppercase()).collect();

    let margin = 8usize;
    let window = adapter.len() + margin;
    let search_start = read_u.len().saturating_sub(window);

    let mut best: Option<AdapterMatch> = None;
    let k_start = search_start.saturating_sub(adapter.len());
    for k in k_start..read_u.len() {
        let overlap = (read_u.len() - k).min(adapter.len());
        if overlap < min_overlap {
            continue;
        }
        if let Some(m) = align_placement_generic(&read_u, &adapter, k, max_error_rate, min_overlap)
        {
            best = Some(match best {
                None => m,
                Some(prev) => prefer(prev, m),
            });
        }
    }

    best
}

/// Align adapter starting at `k` on the read using generic heap-allocated DP.
pub fn align_placement_generic(
    read: &[u8],
    adapter: &[u8],
    k: usize,
    max_error_rate: f32,
    min_overlap: usize,
) -> Option<AdapterMatch> {
    let max_i = (read.len() - k).min(adapter.len());
    if max_i < min_overlap {
        return None;
    }

    let mut mismatches = 0usize;
    for i in 0..max_i {
        if read[k + i] != adapter[i] {
            mismatches += 1;
        }
    }
    let allowed = (max_i as f32 * max_error_rate).floor() as usize;
    if mismatches <= allowed {
        return Some(AdapterMatch {
            read_start: k,
            read_end: k + max_i,
            adapter_start: 0,
            adapter_end: max_i,
            errors: mismatches,
            overlap: max_i,
        });
    }

    let n = max_i;
    let m = max_i;
    let r = &read[k..k + n];
    let a = &adapter[..m];

    let mut prev: Vec<u32> = (0..=m as u32).collect();
    let mut best_end_j = 0usize;
    let mut best_dist = u32::MAX;

    for i in 1..=n {
        let mut curr = vec![0u32; m + 1];
        curr[0] = i as u32;
        for j in 1..=m {
            let cost = if r[i - 1] == a[j - 1] { 0 } else { 1 };
            curr[j] = (prev[j - 1] + cost).min(prev[j] + 1).min(curr[j - 1] + 1);
        }
        if i == n {
            for (j, &dist) in curr.iter().enumerate().skip(min_overlap) {
                if dist < best_dist {
                    best_dist = dist;
                    best_end_j = j;
                }
            }
        }
        prev = curr;
    }

    if best_dist == u32::MAX {
        return None;
    }
    let overlap = best_end_j;
    let errors = best_dist as usize;
    let allowed = (overlap as f32 * max_error_rate).floor() as usize;
    if errors <= allowed && overlap >= min_overlap {
        Some(AdapterMatch {
            read_start: k,
            read_end: k + n,
            adapter_start: 0,
            adapter_end: best_end_j,
            errors,
            overlap,
        })
    } else {
        None
    }
}

/// Zero-heap-allocation banded semi-global alignment using stack arrays.
#[inline]
pub fn align_placement_stack(
    read: &[u8],
    adapter: &[u8],
    k: usize,
    max_error_rate: f32,
    min_overlap: usize,
) -> Option<AdapterMatch> {
    let max_i = (read.len() - k).min(adapter.len());
    if max_i < min_overlap {
        return None;
    }

    // Ungapped check first with early reject.
    let allowed = (max_i as f32 * max_error_rate).floor() as usize;
    let mut mismatches = 0usize;
    let mut ungapped_ok = true;
    for i in 0..max_i {
        if !read[k + i].eq_ignore_ascii_case(&adapter[i]) {
            mismatches += 1;
            if mismatches > allowed {
                ungapped_ok = false;
                break;
            }
        }
    }
    if ungapped_ok {
        return Some(AdapterMatch {
            read_start: k,
            read_end: k + max_i,
            adapter_start: 0,
            adapter_end: max_i,
            errors: mismatches,
            overlap: max_i,
        });
    }

    // If longer than stack buffer, fallback to generic
    if max_i > 63 {
        return align_placement_generic(read, adapter, k, max_error_rate, min_overlap);
    }

    let n = max_i;
    let m = max_i;
    let r = &read[k..k + n];
    let a = &adapter[..m];

    // Ukkonen-style abort when the whole row exceeds any acceptable end distance.
    let dp_abort = ((m as f32 * max_error_rate).floor() as u32).saturating_add(1);

    let mut prev = [0u32; 64];
    for (j, val) in prev.iter_mut().enumerate().take(m + 1) {
        *val = j as u32;
    }

    let mut best_end_j = 0usize;
    let mut best_dist = u32::MAX;

    let mut curr = [0u32; 64];
    for i in 1..=n {
        curr[0] = i as u32;
        let r_byte = r[i - 1].to_ascii_uppercase();
        let mut row_min = curr[0];
        for j in 1..=m {
            let cost = if r_byte == a[j - 1] { 0 } else { 1 };
            let v = (prev[j - 1] + cost).min(prev[j] + 1).min(curr[j - 1] + 1);
            curr[j] = v;
            if v < row_min {
                row_min = v;
            }
        }
        if row_min > dp_abort {
            return None;
        }
        if i == n {
            for (j, &dist) in curr[..=m].iter().enumerate().skip(min_overlap) {
                if dist < best_dist {
                    best_dist = dist;
                    best_end_j = j;
                }
            }
        }
        prev[..=m].copy_from_slice(&curr[..=m]);
    }

    if best_dist == u32::MAX {
        return None;
    }
    let overlap = best_end_j;
    let errors = best_dist as usize;
    let allowed = (overlap as f32 * max_error_rate).floor() as usize;
    if errors <= allowed && overlap >= min_overlap {
        Some(AdapterMatch {
            read_start: k,
            read_end: k + n,
            adapter_start: 0,
            adapter_end: best_end_j,
            errors,
            overlap,
        })
    } else {
        None
    }
}

#[inline]
fn push_candidate(
    cands: &mut [usize; 64],
    n: &mut usize,
    k: usize,
    k_start: usize,
    read_len: usize,
) -> bool {
    if k >= k_start && k < read_len && !cands[..*n].contains(&k) {
        if *n < 64 {
            cands[*n] = k;
            *n += 1;
            true
        } else {
            false
        }
    } else {
        true
    }
}

/// Fast seed-and-extend zero-allocation 3' adapter search.
pub fn find_3p_adapter_fast(read: &[u8], adapter: &CompiledAdapter) -> Option<AdapterMatch> {
    if adapter.seq.is_empty() || read.len() < adapter.min_overlap || adapter.min_overlap == 0 {
        return None;
    }

    // Safety fallback for abnormally long adapters (> 63 bp)
    if adapter.seq.len() > 63 {
        return find_3p_adapter_generic(
            read,
            &adapter.seq,
            adapter.max_error_rate,
            adapter.min_overlap,
        );
    }

    let margin = 8usize;
    let window = adapter.seq.len() + margin;
    let search_start = read.len().saturating_sub(window);
    let k_start = search_start.saturating_sub(adapter.seq.len());

    let mut best: Option<AdapterMatch> = None;

    // Fast Path 1: Ungapped scan across all placements with early break
    for k in k_start..read.len() {
        let overlap = (read.len() - k).min(adapter.seq.len());
        if overlap < adapter.min_overlap {
            continue;
        }
        let allowed = (overlap as f32 * adapter.max_error_rate).floor() as usize;
        let mut mismatches = 0usize;
        let mut exceeded = false;
        for i in 0..overlap {
            if !read[k + i].eq_ignore_ascii_case(&adapter.seq[i]) {
                mismatches += 1;
                if mismatches > allowed {
                    exceeded = true;
                    break;
                }
            }
        }
        if !exceeded {
            let m = AdapterMatch {
                read_start: k,
                read_end: k + overlap,
                adapter_start: 0,
                adapter_end: overlap,
                errors: mismatches,
                overlap,
            };
            best = Some(match best {
                None => m,
                Some(prev) => prefer(prev, m),
            });
        }
    }

    // Early exit if an exact full match was found
    if let Some(ref m) = best {
        if m.errors == 0 && m.overlap == adapter.seq.len() {
            return best;
        }
    }

    // Fast Path 2: inverted 8-mer seed lookup — scan read tail once, O(1) probe.
    let mut num_candidates = 0usize;
    let mut candidate_positions = [0usize; 64];
    let mut overflowed = false;

    if read.len() >= 8 && !adapter.seed_lookup.is_empty() {
        let max_k = read.len() - 8;
        let mut hit_offsets = [0u8; 4];

        for p in k_start..=max_k {
            let chunk = &read[p..p + 8];
            let word = u64::from_ne_bytes(chunk.try_into().unwrap()) & !0x2020_2020_2020_2020;
            let n_hits = adapter.seed_lookup.lookup(word, &mut hit_offsets);
            for &off in hit_offsets.iter().take(n_hits) {
                let offset = off as usize;
                let center = p.saturating_sub(offset);
                for delta in -4isize..=4 {
                    let k = center as isize + delta;
                    if k < 0 {
                        continue;
                    }
                    let k = k as usize;
                    if !push_candidate(
                        &mut candidate_positions,
                        &mut num_candidates,
                        k,
                        k_start,
                        read.len(),
                    ) {
                        overflowed = true;
                        break;
                    }
                }
                if overflowed {
                    break;
                }
            }
            if overflowed {
                break;
            }
        }

        if overflowed {
            return find_3p_adapter_generic(
                read,
                &adapter.seq,
                adapter.max_error_rate,
                adapter.min_overlap,
            );
        }

        for &k in &candidate_positions[..num_candidates] {
            if let Some(m) = align_placement_stack(
                read,
                &adapter.seq,
                k,
                adapter.max_error_rate,
                adapter.min_overlap,
            ) {
                best = Some(match best {
                    None => m,
                    Some(prev) => prefer(prev, m),
                });
            }
        }
    }

    // Fast Path 3: short-overlap indel rescue.
    // Unconditional DP over the full 31-base tail was the Sprint 2.5→3 regression
    // (~29 stack-DP verifications per negative read). Keep correctness for short
    // indel windows via 4-mer candidate filtering; only tiny overlaps (<4) are
    // scanned exhaustively.
    if adapter.seq.len() < 16 {
        for k in k_start..read.len() {
            let overlap = (read.len() - k).min(adapter.seq.len());
            if overlap < adapter.min_overlap {
                continue;
            }
            if let Some(m) = align_placement_stack(
                read,
                &adapter.seq,
                k,
                adapter.max_error_rate,
                adapter.min_overlap,
            ) {
                best = Some(match best {
                    None => m,
                    Some(prev) => prefer(prev, m),
                });
            }
        }
    } else {
        // Tiny overlaps cannot host a 4-mer seed.
        let tiny_start = read.len().saturating_sub(3).max(k_start);
        for k in tiny_start..read.len() {
            let overlap = read.len() - k;
            if overlap < adapter.min_overlap {
                continue;
            }
            if candidate_positions[..num_candidates].contains(&k) {
                continue;
            }
            if let Some(m) = align_placement_stack(
                read,
                &adapter.seq,
                k,
                adapter.max_error_rate,
                adapter.min_overlap,
            ) {
                best = Some(match best {
                    None => m,
                    Some(prev) => prefer(prev, m),
                });
            }
        }

        // 4-mer inverted scan of the short 3' window (overlap < 16).
        if adapter.short4_len > 0 && read.len() >= 4 {
            let win_start = read.len().saturating_sub(15).max(k_start);
            let max_p = read.len() - 4;
            if win_start <= max_p {
                let mut short_cands = [0usize; 64];
                let mut n_short = 0usize;
                for p in win_start..=max_p {
                    let word =
                        u32::from_ne_bytes(read[p..p + 4].try_into().unwrap()) & !0x2020_2020;
                    for i in 0..adapter.short4_len as usize {
                        let (seed, offset) = adapter.short4[i];
                        if word != seed {
                            continue;
                        }
                        let center = p.saturating_sub(offset as usize);
                        for delta in -2isize..=2 {
                            let k = center as isize + delta;
                            if k < 0 {
                                continue;
                            }
                            let k = k as usize;
                            if k < win_start || k >= read.len() {
                                continue;
                            }
                            if candidate_positions[..num_candidates].contains(&k) {
                                continue;
                            }
                            let _ = push_candidate(
                                &mut short_cands,
                                &mut n_short,
                                k,
                                win_start,
                                read.len(),
                            );
                        }
                    }
                }
                for &k in &short_cands[..n_short] {
                    if let Some(m) = align_placement_stack(
                        read,
                        &adapter.seq,
                        k,
                        adapter.max_error_rate,
                        adapter.min_overlap,
                    ) {
                        best = Some(match best {
                            None => m,
                            Some(prev) => prefer(prev, m),
                        });
                    }
                }
            }
        }
    }

    best
}

/// Find best 3' adapter match in the read.
/// Precompiles adapter on the fly and invokes the fast path.
pub fn find_3p_adapter(
    read: &[u8],
    adapter: &[u8],
    max_error_rate: f32,
    min_overlap: usize,
) -> Option<AdapterMatch> {
    if adapter.is_empty() || read.is_empty() || min_overlap == 0 {
        return None;
    }
    let compiled = CompiledAdapter::new(adapter, max_error_rate, min_overlap);
    compiled.find_3p(read)
}

/// Trim 3' adapter from read if found; returns (new_seq_end, was_trimmed, overlap).
pub fn trim_3p_adapter(
    seq: &[u8],
    adapter: &[u8],
    max_error_rate: f32,
    min_overlap: usize,
) -> (usize, bool, usize) {
    if adapter.is_empty() {
        return (seq.len(), false, 0);
    }
    let compiled = CompiledAdapter::new(adapter, max_error_rate, min_overlap);
    compiled.trim_3p(seq)
}

#[cfg(test)]
mod tests {
    use super::*;

    const ILLUMINA_R1: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
    const ILLUMINA_R2: &[u8] = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

    #[test]
    fn exact_full_adapter() {
        let read = b"ACGTACGTAGATCGGAAG";
        let adapter = b"AGATCGGAAG";
        let m = find_3p_adapter(read, adapter, 0.1, 3).unwrap();
        assert_eq!(m.read_start, 8);
        assert_eq!(m.errors, 0);
    }

    #[test]
    fn partial_adapter() {
        let read = b"ACGTACGTAGATC";
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
        let m = find_3p_adapter(read, adapter, 0.1, 3).unwrap();
        assert!(m.overlap >= 3);
        assert_eq!(&read[m.read_start..], b"AGATC");
    }

    #[test]
    fn no_adapter() {
        let read = b"AAAAAAAAAAAAAAAA";
        let adapter = b"AGATCGGAAGAGC";
        assert!(find_3p_adapter(read, adapter, 0.1, 5).is_none());
    }

    #[test]
    fn trim_helper() {
        let seq = b"ACGTAGATCGG";
        let (end, trimmed, _) = trim_3p_adapter(seq, b"AGATCGG", 0.1, 3);
        assert!(trimmed);
        assert_eq!(end, 4);
    }

    // =========================================================================
    // Differential Tests: fast adapter path vs existing generic path
    // =========================================================================

    fn assert_fast_eq_generic(
        read: &[u8],
        adapter: &[u8],
        max_error_rate: f32,
        min_overlap: usize,
    ) {
        let generic_res = find_3p_adapter_generic(read, adapter, max_error_rate, min_overlap);
        let compiled = CompiledAdapter::new(adapter, max_error_rate, min_overlap);
        let fast_res = find_3p_adapter_fast(read, &compiled);
        assert_eq!(
            fast_res,
            generic_res,
            "Differential mismatch!\nRead: {}\nAdapter: {}\nFast: {:?}\nGeneric: {:?}",
            String::from_utf8_lossy(read),
            String::from_utf8_lossy(adapter),
            fast_res,
            generic_res
        );
    }

    #[test]
    fn diff_full_exact_adapter() {
        let prefixes = [
            b"ACGTACGT".as_slice(),
            b"TTTTGGGGCCCCAAAA".as_slice(),
            b"".as_slice(),
            b"NNNNNNNNNNNNNNNN".as_slice(),
        ];
        for prefix in prefixes {
            let mut read = prefix.to_vec();
            read.extend_from_slice(ILLUMINA_R1);
            assert_fast_eq_generic(&read, ILLUMINA_R1, 0.1, 3);
            assert_fast_eq_generic(&read, ILLUMINA_R1, 0.2, 5);

            let mut read2 = prefix.to_vec();
            read2.extend_from_slice(ILLUMINA_R2);
            assert_fast_eq_generic(&read2, ILLUMINA_R2, 0.1, 3);
        }
    }

    #[test]
    fn diff_partial_3p_adapter() {
        let prefix = b"ACGTACGTACGTACGT";
        for overlap in 1..=ILLUMINA_R1.len() {
            let mut read = prefix.to_vec();
            read.extend_from_slice(&ILLUMINA_R1[..overlap]);
            for min_overlap in [3, 5, 8, 10, 15] {
                assert_fast_eq_generic(&read, ILLUMINA_R1, 0.1, min_overlap);
            }
        }
    }

    #[test]
    fn diff_no_adapter() {
        let reads = [
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".as_slice(),
            b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".as_slice(),
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".as_slice(),
            b"TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA".as_slice(),
            b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".as_slice(),
        ];
        for read in reads {
            assert_fast_eq_generic(read, ILLUMINA_R1, 0.1, 3);
            assert_fast_eq_generic(read, ILLUMINA_R2, 0.1, 3);
        }
    }

    #[test]
    fn diff_adapter_with_allowed_mismatch() {
        let prefix = b"ACGTACGTACGTACGT";
        // 34 bp adapter with error rate 0.1 allows floor(34 * 0.1) = 3 mismatches
        let mut ad_1mis = ILLUMINA_R1.to_vec();
        ad_1mis[2] = b'T'; // A -> T
        let mut read = prefix.to_vec();
        read.extend_from_slice(&ad_1mis);
        assert_fast_eq_generic(&read, ILLUMINA_R1, 0.1, 3);

        let mut ad_2mis = ILLUMINA_R1.to_vec();
        ad_2mis[2] = b'T';
        ad_2mis[10] = b'A';
        let mut read2 = prefix.to_vec();
        read2.extend_from_slice(&ad_2mis);
        assert_fast_eq_generic(&read2, ILLUMINA_R1, 0.1, 3);

        let mut ad_3mis = ILLUMINA_R1.to_vec();
        ad_3mis[2] = b'T';
        ad_3mis[10] = b'A';
        ad_3mis[20] = b'C';
        let mut read3 = prefix.to_vec();
        read3.extend_from_slice(&ad_3mis);
        assert_fast_eq_generic(&read3, ILLUMINA_R1, 0.1, 3);
    }

    #[test]
    fn diff_adapter_outside_allowed_error_threshold() {
        let prefix = b"ACGTACGTACGTACGT";
        // 4 mismatches > allowed 3
        let mut ad_4mis = ILLUMINA_R1.to_vec();
        ad_4mis[2] = b'T';
        ad_4mis[10] = b'A';
        ad_4mis[20] = b'C';
        ad_4mis[30] = b'G';
        let mut read = prefix.to_vec();
        read.extend_from_slice(&ad_4mis);
        assert_fast_eq_generic(&read, ILLUMINA_R1, 0.1, 3);
    }

    #[test]
    fn diff_short_overlap_boundary_cases() {
        let prefix = b"ACGTACGTACGT";
        // min_overlap = 5
        for len in 1..=8 {
            let mut read = prefix.to_vec();
            read.extend_from_slice(&ILLUMINA_R1[..len]);
            assert_fast_eq_generic(&read, ILLUMINA_R1, 0.1, 5);
        }
        // min_overlap = 3
        for len in 1..=5 {
            let mut read = prefix.to_vec();
            read.extend_from_slice(&ILLUMINA_R1[..len]);
            assert_fast_eq_generic(&read, ILLUMINA_R1, 0.1, 3);
        }
    }

    #[test]
    fn diff_adapter_with_indels() {
        let prefix = b"ACGTACGTACGT";
        // Insertion in read (adapter deletion)
        let mut ad_del = ILLUMINA_R1[..10].to_vec();
        ad_del.extend_from_slice(&ILLUMINA_R1[11..]); // omit base 10
        let mut read = prefix.to_vec();
        read.extend_from_slice(&ad_del);
        assert_fast_eq_generic(&read, ILLUMINA_R1, 0.1, 3);

        // Deletion in read (adapter insertion)
        let mut ad_ins = ILLUMINA_R1[..10].to_vec();
        ad_ins.push(b'N'); // inserted base
        ad_ins.extend_from_slice(&ILLUMINA_R1[10..]);
        let mut read2 = prefix.to_vec();
        read2.extend_from_slice(&ad_ins);
        assert_fast_eq_generic(&read2, ILLUMINA_R1, 0.1, 3);
    }

    #[test]
    fn diff_short_adapters() {
        let short_ad = b"AGATCGG";
        let prefix = b"ACGTACGTACGT";
        for overlap in 1..=short_ad.len() {
            let mut read = prefix.to_vec();
            read.extend_from_slice(&short_ad[..overlap]);
            assert_fast_eq_generic(&read, short_ad, 0.1, 3);
        }
    }
}
