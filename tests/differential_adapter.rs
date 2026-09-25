use seqmux::trim::{find_3p_adapter_fast, find_3p_adapter_generic, CompiledAdapter};

const ILLUMINA_R1: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
const ILLUMINA_R2: &[u8] = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

fn assert_differential(
    read: &[u8],
    adapter: &[u8],
    max_error_rate: f32,
    min_overlap: usize,
    context: &str,
) {
    let generic = find_3p_adapter_generic(read, adapter, max_error_rate, min_overlap);
    let compiled = CompiledAdapter::new(adapter, max_error_rate, min_overlap);
    let fast = find_3p_adapter_fast(read, &compiled);
    assert_eq!(
        fast,
        generic,
        "Differential mismatch [{}]:\nRead:    {}\nAdapter: {}\nFast:    {:?}\nGeneric: {:?}",
        context,
        String::from_utf8_lossy(read),
        String::from_utf8_lossy(adapter),
        fast,
        generic,
    );
}

#[test]
fn test_differential_full_exact_adapter() {
    let prefixes: [&[u8]; 6] = [
        b"",
        b"ACGT",
        b"ACGTACGT",
        b"TTTTGGGGCCCCAAAA",
        b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        b"NNNNNNNNNNNNNNNN",
    ];

    for prefix in prefixes {
        let mut r1 = prefix.to_vec();
        r1.extend_from_slice(ILLUMINA_R1);
        assert_differential(&r1, ILLUMINA_R1, 0.1, 3, "full exact R1");
        assert_differential(&r1, ILLUMINA_R1, 0.2, 5, "full exact R1 error=0.2");

        let mut r2 = prefix.to_vec();
        r2.extend_from_slice(ILLUMINA_R2);
        assert_differential(&r2, ILLUMINA_R2, 0.1, 3, "full exact R2");
    }
}

#[test]
fn test_differential_partial_3p_adapter() {
    let prefix = b"ACGTACGTACGTACGT";
    for overlap in 1..=ILLUMINA_R1.len() {
        let mut read = prefix.to_vec();
        read.extend_from_slice(&ILLUMINA_R1[..overlap]);
        for min_overlap in [1, 2, 3, 5, 8, 10, 15] {
            assert_differential(
                &read,
                ILLUMINA_R1,
                0.1,
                min_overlap,
                &format!("partial overlap={overlap} min_overlap={min_overlap}"),
            );
        }
    }
}

#[test]
fn test_differential_no_adapter() {
    let reads: [&[u8]; 6] = [
        b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",
        b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
    ];

    for read in reads {
        assert_differential(read, ILLUMINA_R1, 0.1, 3, "no adapter R1");
        assert_differential(read, ILLUMINA_R2, 0.1, 3, "no adapter R2");
    }
}

#[test]
fn test_differential_adapter_with_allowed_mismatch() {
    let prefix = b"ACGTACGTACGTACGT";

    // 1 mismatch
    for pos in [0, 5, 12, 20, 33] {
        let mut ad = ILLUMINA_R1.to_vec();
        ad[pos] = match ad[pos] {
            b'A' => b'T',
            _ => b'A',
        };
        let mut read = prefix.to_vec();
        read.extend_from_slice(&ad);
        assert_differential(&read, ILLUMINA_R1, 0.1, 3, &format!("1 mismatch at {pos}"));
    }

    // 2 mismatches
    let mut ad2 = ILLUMINA_R1.to_vec();
    ad2[2] = b'T';
    ad2[15] = b'G';
    let mut read2 = prefix.to_vec();
    read2.extend_from_slice(&ad2);
    assert_differential(&read2, ILLUMINA_R1, 0.1, 3, "2 mismatches (2, 15)");

    // 3 mismatches (allowed: floor(34 * 0.1) = 3)
    let mut ad3 = ILLUMINA_R1.to_vec();
    ad3[2] = b'T';
    ad3[15] = b'G';
    ad3[28] = b'C';
    let mut read3 = prefix.to_vec();
    read3.extend_from_slice(&ad3);
    assert_differential(&read3, ILLUMINA_R1, 0.1, 3, "3 mismatches (2, 15, 28)");
}

#[test]
fn test_differential_adapter_outside_allowed_error_threshold() {
    let prefix = b"ACGTACGTACGTACGT";

    // 4 mismatches (threshold is 3 for len 34, error 0.1)
    let mut ad4 = ILLUMINA_R1.to_vec();
    ad4[2] = b'T';
    ad4[10] = b'G';
    ad4[20] = b'C';
    ad4[30] = b'A';
    let mut read = prefix.to_vec();
    read.extend_from_slice(&ad4);
    assert_differential(
        &read,
        ILLUMINA_R1,
        0.1,
        3,
        "4 mismatches (exceeds threshold)",
    );

    // Partial overlap with error exceeding threshold: overlap 10, error rate 0.1 allows 1 mismatch
    let mut ad_part = ILLUMINA_R1[..10].to_vec();
    ad_part[2] = b'T';
    ad_part[8] = b'G'; // 2 mismatches in 10 bp > floor(1)
    let mut read_part = prefix.to_vec();
    read_part.extend_from_slice(&ad_part);
    assert_differential(
        &read_part,
        ILLUMINA_R1,
        0.1,
        3,
        "partial 10 bp with 2 mismatches",
    );
}

#[test]
fn test_differential_short_overlap_boundary_cases() {
    let prefix = b"ACGTACGTACGTACGT";
    for min_overlap in [1, 2, 3, 4, 5, 8] {
        for len in 0..=10 {
            let mut read = prefix.to_vec();
            read.extend_from_slice(&ILLUMINA_R1[..len]);
            assert_differential(
                &read,
                ILLUMINA_R1,
                0.1,
                min_overlap,
                &format!("boundary len={len} min_overlap={min_overlap}"),
            );
        }
    }
}

#[test]
fn test_differential_indels() {
    let prefix = b"ACGTACGTACGTACGT";

    // Single deletion in adapter (insertion in read)
    for del_pos in [2, 7, 14, 22] {
        let mut ad_del = ILLUMINA_R1[..del_pos].to_vec();
        ad_del.extend_from_slice(&ILLUMINA_R1[del_pos + 1..]);
        let mut read = prefix.to_vec();
        read.extend_from_slice(&ad_del);
        assert_differential(
            &read,
            ILLUMINA_R1,
            0.1,
            3,
            &format!("deletion at adapter pos {del_pos}"),
        );
    }

    // Single insertion in adapter (deletion in read)
    for ins_pos in [2, 7, 14, 22] {
        let mut ad_ins = ILLUMINA_R1[..ins_pos].to_vec();
        ad_ins.push(b'N');
        ad_ins.extend_from_slice(&ILLUMINA_R1[ins_pos..]);
        let mut read = prefix.to_vec();
        read.extend_from_slice(&ad_ins);
        assert_differential(
            &read,
            ILLUMINA_R1,
            0.1,
            3,
            &format!("insertion at adapter pos {ins_pos}"),
        );
    }
}

#[test]
fn test_differential_synthetic_random_stress() {
    // Deterministic pseudo-random sequence generator (LCG)
    let mut state = 123456789u64;
    let mut lcg = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        state
    };
    let bases = *b"ACGT";

    for i in 0..5000 {
        let read_len = 50 + (lcg() % 50) as usize;
        let mut read: Vec<u8> = (0..read_len).map(|_| bases[(lcg() % 4) as usize]).collect();

        // 30% chance to append a modified adapter
        if lcg() % 100 < 30 {
            let overlap = 3 + (lcg() % (ILLUMINA_R1.len() - 2) as u64) as usize;
            let mut ad = ILLUMINA_R1[..overlap].to_vec();
            // maybe mutate a base
            if lcg() % 100 < 50 && !ad.is_empty() {
                let m_pos = (lcg() as usize) % ad.len();
                ad[m_pos] = bases[(lcg() % 4) as usize];
            }
            read.extend_from_slice(&ad);
        }

        assert_differential(
            &read,
            ILLUMINA_R1,
            0.1,
            3,
            &format!("random synthetic read #{i}"),
        );
    }
}
