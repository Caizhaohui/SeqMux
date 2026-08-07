/// 3' adapter match result.
#[derive(Debug, Clone, PartialEq)]
pub struct AdapterMatch {
    pub read_start: usize,
    pub read_end: usize,
    pub adapter_start: usize,
    pub adapter_end: usize,
    pub errors: usize,
    pub overlap: usize,
}

/// Find best 3' adapter match in the read.
///
/// Algorithm:
/// 1. Prefer exact suffix/prefix overlaps.
/// 2. Run small semi-global DP over the read tail allowing substitutions and indels.
/// 3. Accept if errors <= floor(overlap * max_error_rate) and overlap >= min_overlap.
pub fn find_3p_adapter(
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

    // Search window: adapter length + margin
    let margin = 8usize;
    let window = adapter.len() + margin;
    let search_start = read_u.len().saturating_sub(window);

    let mut best: Option<AdapterMatch> = None;

    // Evaluate all placements where adapter[0..] aligns starting at read position k
    // (semi-global: adapter can hang off the 3' end of the read).
    // k can be from search_start to read_u.len()-1 (partial suffix of adapter).
    // Also allow adapter starting before search_start if full adapter fits near end.
    let k_start = search_start.saturating_sub(adapter.len());
    for k in k_start..read_u.len() {
        // Align adapter starting at k; read may end before adapter ends.
        let overlap = (read_u.len() - k).min(adapter.len());
        if overlap < min_overlap {
            continue;
        }
        // DP for this placement with free end gaps only at adapter 3' and read 5' of window
        if let Some(m) = align_placement(&read_u, &adapter, k, max_error_rate, min_overlap) {
            best = Some(match best {
                None => m,
                Some(prev) => prefer(prev, m),
            });
        }
    }

    best
}

fn prefer(a: AdapterMatch, b: AdapterMatch) -> AdapterMatch {
    // Prefer more overlap, then fewer errors, then earlier start
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

/// Align adapter starting at `k` on the read using banded edit distance (NW-style).
fn align_placement(
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

    // Simple approach: count mismatches for ungapped, then try small indel DP
    // Ungapped first (fast path)
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

    // Semi-global DP for indels when ungapped fails
    // Align read[k..] with adapter[0..], free gaps at end of shorter side for partial
    let n = max_i; // use only overlapping portion length upper bound
    let m = max_i; // adapter prefix of same length
    let r = &read[k..k + n];
    let a = &adapter[..m];

    // dp[i][j] = edit distance for r[0..i) vs a[0..j)
    let mut prev: Vec<u32> = (0..=m as u32).collect();
    // Actually for 3' adapter: free end gaps on adapter only when partial at read end.
    // Standard cutadapt semi-global: no penalty for adapter bases beyond read end
    // (already handled by limiting m). Start gaps on read free? No for 3'.

    let mut best_end_j = 0usize;
    let mut best_dist = u32::MAX;

    for i in 1..=n {
        let mut curr = vec![0u32; m + 1];
        curr[0] = i as u32;
        for j in 1..=m {
            let cost = if r[i - 1] == a[j - 1] { 0 } else { 1 };
            curr[j] = (prev[j - 1] + cost).min(prev[j] + 1).min(curr[j - 1] + 1);
        }
        // At end of read, take best over adapter prefixes with enough overlap
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
            read_end: k + n, // trim from k to end of read for 3' adapter
            adapter_start: 0,
            adapter_end: best_end_j,
            errors,
            overlap,
        })
    } else {
        None
    }
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
    match find_3p_adapter(seq, adapter, max_error_rate, min_overlap) {
        Some(m) => (m.read_start, true, m.overlap),
        None => (seq.len(), false, 0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
