use seqmux::trim::{find_3p_adapter_fast, find_3p_adapter_generic, CompiledAdapter};

const ILLUMINA_R1: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
#[allow(dead_code)]
const ILLUMINA_R2: &[u8] = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
const BASES: [u8; 4] = *b"ACGT";

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
        fast, generic,
        "Adversarial divergence [{context}]:\nRead:    {}\nAdapter: {}\nFast:    {:?}\nGeneric: {:?}",
        String::from_utf8_lossy(read),
        String::from_utf8_lossy(adapter),
        fast, generic
    );
}

#[test]
fn test_exhaustive_all_single_substitutions() {
    let prefix = b"ACGTACGTACGTACGT";
    // Test every single substitution position on ILLUMINA_R1 (34 bp)
    for pos in 0..ILLUMINA_R1.len() {
        for &sub in &BASES {
            if sub == ILLUMINA_R1[pos] {
                continue;
            }
            let mut ad = ILLUMINA_R1.to_vec();
            ad[pos] = sub;

            let mut read = prefix.to_vec();
            read.extend_from_slice(&ad);
            assert_differential(
                &read,
                ILLUMINA_R1,
                0.1,
                3,
                &format!("1 sub pos={pos} base={}", sub as char),
            );
        }
    }
}

#[test]
fn test_exhaustive_all_single_deletions() {
    let prefix = b"ACGTACGTACGTACGT";
    // Delete base at every single position on ILLUMINA_R1
    for del_pos in 0..ILLUMINA_R1.len() {
        let mut ad = ILLUMINA_R1[..del_pos].to_vec();
        ad.extend_from_slice(&ILLUMINA_R1[del_pos + 1..]);

        let mut read = prefix.to_vec();
        read.extend_from_slice(&ad);
        assert_differential(&read, ILLUMINA_R1, 0.1, 3, &format!("1 del pos={del_pos}"));
    }
}

#[test]
fn test_exhaustive_all_single_insertions() {
    let prefix = b"ACGTACGTACGTACGT";
    // Insert base at every single position 0..=34 on ILLUMINA_R1
    for ins_pos in 0..=ILLUMINA_R1.len() {
        for &ins_base in &BASES {
            let mut ad = ILLUMINA_R1[..ins_pos].to_vec();
            ad.push(ins_base);
            ad.extend_from_slice(&ILLUMINA_R1[ins_pos..]);

            let mut read = prefix.to_vec();
            read.extend_from_slice(&ad);
            assert_differential(
                &read,
                ILLUMINA_R1,
                0.1,
                3,
                &format!("1 ins pos={ins_pos} base={}", ins_base as char),
            );
        }
    }
}

#[test]
fn test_exhaustive_two_edit_combinations() {
    let prefix = b"ACGTACGTACGTACGT";

    // 1. Two substitutions: sample across pairs (step by 3 to cover all regions)
    for i in (0..ILLUMINA_R1.len()).step_by(2) {
        for j in (i + 1..ILLUMINA_R1.len()).step_by(2) {
            let mut ad = ILLUMINA_R1.to_vec();
            ad[i] = if ad[i] == b'A' { b'T' } else { b'A' };
            ad[j] = if ad[j] == b'C' { b'G' } else { b'C' };

            let mut read = prefix.to_vec();
            read.extend_from_slice(&ad);
            assert_differential(&read, ILLUMINA_R1, 0.1, 3, &format!("2 subs pos=({i},{j})"));
        }
    }

    // 2. Two deletions: sample across pairs
    for i in (0..ILLUMINA_R1.len()).step_by(3) {
        for j in (i + 1..ILLUMINA_R1.len()).step_by(3) {
            let mut ad = Vec::new();
            for (idx, &b) in ILLUMINA_R1.iter().enumerate() {
                if idx != i && idx != j {
                    ad.push(b);
                }
            }
            let mut read = prefix.to_vec();
            read.extend_from_slice(&ad);
            assert_differential(&read, ILLUMINA_R1, 0.1, 3, &format!("2 dels pos=({i},{j})"));
        }
    }

    // 3. One deletion + one substitution
    for d in (0..ILLUMINA_R1.len()).step_by(3) {
        for s in (0..ILLUMINA_R1.len()).step_by(3) {
            if d == s {
                continue;
            }
            let mut ad = ILLUMINA_R1.to_vec();
            ad[s] = if ad[s] == b'A' { b'T' } else { b'A' };
            ad.remove(d);

            let mut read = prefix.to_vec();
            read.extend_from_slice(&ad);
            assert_differential(
                &read,
                ILLUMINA_R1,
                0.1,
                3,
                &format!("1 del pos={d}, 1 sub pos={s}"),
            );
        }
    }

    // 4. One insertion + one substitution
    for ins in (0..=ILLUMINA_R1.len()).step_by(3) {
        for s in (0..ILLUMINA_R1.len()).step_by(3) {
            let mut ad = ILLUMINA_R1.to_vec();
            ad[s] = if ad[s] == b'A' { b'T' } else { b'A' };
            ad.insert(ins, b'N');

            let mut read = prefix.to_vec();
            read.extend_from_slice(&ad);
            assert_differential(
                &read,
                ILLUMINA_R1,
                0.1,
                3,
                &format!("1 ins pos={ins}, 1 sub pos={s}"),
            );
        }
    }
}

#[test]
fn test_exhaustive_three_edit_combinations() {
    let prefix = b"ACGTACGTACGT";

    // 3 substitutions strategically spaced across the 34 bp adapter to challenge 8-mer seeds:
    // E.g. errors at positions 7, 15, 23 (breaks pieces into 7, 7, 7, 10)
    for p1 in [3, 5, 7, 9] {
        for p2 in [12, 15, 18] {
            for p3 in [21, 24, 27] {
                let mut ad = ILLUMINA_R1.to_vec();
                ad[p1] = if ad[p1] == b'A' { b'C' } else { b'A' };
                ad[p2] = if ad[p2] == b'A' { b'G' } else { b'A' };
                ad[p3] = if ad[p3] == b'A' { b'T' } else { b'A' };

                let mut read = prefix.to_vec();
                read.extend_from_slice(&ad);
                assert_differential(
                    &read,
                    ILLUMINA_R1,
                    0.1,
                    3,
                    &format!("3 subs ({p1},{p2},{p3})"),
                );
            }
        }
    }

    // Mixed 3 edits: 1 deletion + 2 substitutions
    for del in [4, 12, 20] {
        for s1 in [2, 10, 18] {
            for s2 in [16, 24, 30] {
                let mut ad = ILLUMINA_R1.to_vec();
                ad[s1] = b'T';
                ad[s2] = b'G';
                ad.remove(del);

                let mut read = prefix.to_vec();
                read.extend_from_slice(&ad);
                assert_differential(
                    &read,
                    ILLUMINA_R1,
                    0.1,
                    3,
                    &format!("1 del {del}, 2 subs ({s1},{s2})"),
                );
            }
        }
    }
}

#[test]
fn test_exhaustive_partial_overlaps_with_indels() {
    let prefix = b"ACGTACGTACGT";

    // For every partial overlap length from 3 to 34 bp:
    for overlap in 3..=ILLUMINA_R1.len() {
        let partial = &ILLUMINA_R1[..overlap];

        // 1. Exact partial
        let mut read = prefix.to_vec();
        read.extend_from_slice(partial);
        assert_differential(
            &read,
            ILLUMINA_R1,
            0.1,
            3,
            &format!("exact partial overlap={overlap}"),
        );

        // 2. Partial with 1 substitution at every possible position
        for s in 0..overlap {
            let mut part_mut = partial.to_vec();
            part_mut[s] = if part_mut[s] == b'A' { b'T' } else { b'A' };
            let mut read_mut = prefix.to_vec();
            read_mut.extend_from_slice(&part_mut);
            assert_differential(
                &read_mut,
                ILLUMINA_R1,
                0.1,
                3,
                &format!("partial {overlap} sub {s}"),
            );
        }

        // 3. Partial with 1 deletion at every possible position (if overlap > 3)
        if overlap > 3 {
            for d in 0..overlap {
                let mut part_del = partial.to_vec();
                part_del.remove(d);
                let mut read_del = prefix.to_vec();
                read_del.extend_from_slice(&part_del);
                assert_differential(
                    &read_del,
                    ILLUMINA_R1,
                    0.1,
                    3,
                    &format!("partial {overlap} del {d}"),
                );
            }
        }

        // 4. Partial with 1 insertion at every possible position
        for ins in 0..=overlap {
            let mut part_ins = partial.to_vec();
            part_ins.insert(ins, b'N');
            let mut read_ins = prefix.to_vec();
            read_ins.extend_from_slice(&part_ins);
            assert_differential(
                &read_ins,
                ILLUMINA_R1,
                0.1,
                3,
                &format!("partial {overlap} ins {ins}"),
            );
        }
    }
}

#[test]
fn test_all_min_overlap_boundaries() {
    let prefix = b"ACGTACGTACGT";
    let test_overlaps = [1, 2, 3, 4, 5, 8, 10, 15];

    for &min_overlap in &test_overlaps {
        for len in 0..=(min_overlap + 3) {
            if len <= ILLUMINA_R1.len() {
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
}

#[test]
fn test_property_based_random_adversarial_stress() {
    let mut state = 987654321012345u64;
    let mut lcg = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        state
    };

    for iter in 0..10_000 {
        let prefix_len = (lcg() % 50) as usize;
        let mut read: Vec<u8> = (0..prefix_len)
            .map(|_| BASES[(lcg() % 4) as usize])
            .collect();

        // Randomly pick overlap length 3..=34
        let overlap = 3 + (lcg() % (ILLUMINA_R1.len() - 2) as u64) as usize;
        let mut ad = ILLUMINA_R1[..overlap].to_vec();

        // Number of edits: 0, 1, 2, 3, or 4
        let n_edits = (lcg() % 5) as usize;
        for _ in 0..n_edits {
            if ad.is_empty() {
                break;
            }
            let edit_type = lcg() % 3;
            match edit_type {
                0 => {
                    // Substitution
                    let pos = (lcg() as usize) % ad.len();
                    ad[pos] = BASES[(lcg() % 4) as usize];
                }
                1 => {
                    // Deletion
                    let pos = (lcg() as usize) % ad.len();
                    ad.remove(pos);
                }
                _ => {
                    // Insertion
                    let pos = (lcg() as usize) % (ad.len() + 1);
                    ad.insert(pos, BASES[(lcg() % 4) as usize]);
                }
            }
        }

        read.extend_from_slice(&ad);
        assert_differential(
            &read,
            ILLUMINA_R1,
            0.1,
            3,
            &format!("random property stress #{iter}"),
        );
    }
}
