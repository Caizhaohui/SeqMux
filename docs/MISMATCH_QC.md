# Mismatch=1 QC

Default SeqMux mismatch is **0**, matching the lab Python demux.

## Decision

**Keep default `--mismatches-1/2 = 0`.** Dual 8-mer pairs in I464 are well separated (minimum combined Hamming between different samples = 8), so mismatch=1 does not jump to another sample. Extra hits recover substitution errors that Python exact demux leaves unassigned; they remain opt-in.

Use `--mismatches-1 1 --mismatches-2 1` only when the table’s minimum combined Hamming is ≥ 3 (I464 is 8) **and** you accept recovering those extra pairs.

## I464 200k-pair subsample

| | mm=0 | mm=1 |
|--|------|------|
| assigned | 124,867 (62.43%) | 130,644 (65.32%) |
| extra vs exact | — | **+5,777** |
| ambiguous | 0 | 0 |
| chimeric pairs absorbed | — | **0** / 19,027 chimeric exact-piece pairs |

Of the +5,777:

- 5,537 distance 1 (one barcode off by 1)
- 240 distance 2 (each barcode off by 1)
- largest extra: I464469-D7 (+875), matching the known `GGAACGTT` → `GGAACGTA` tail error

SeqMux mm=1 counts match the independent auditor (`scripts/audit_mismatch1.py`).

## I464 full (140,771,720 pairs, 2026-08-13)

| | mm=0 | mm=1 |
|--|------|------|
| assigned | 87,477,095 (62.141%) | 90,802,891 (64.504%) |
| extra | — | **+3,325,796** (2.36% of total) |
| ambiguous | 0 | 0 |
| any sample count decreased | — | **no** |
| wall (16 threads) | 5:58 | 6:06 |

Top sample extras under mm=1:

| Sample | Extra |
|--------|-------|
| I464469-D7 | +661,193 |
| I464469-B1 | +127,450 |
| I464469-C6 | +124,689 |
| I464469-A6 | +121,081 |
| I464469-C4 | +110,922 |

D7 dominates because many unmatched pairs are 1-edit from `CGCTATGT+GGAACGTT` (see Python `unmatched_barcode_counts.csv`).

Barcode design: min BC1–BC1 Hamming = 2, min BC2–BC2 = 3, but **combined sample-pair Hamming ≥ 8**, so a single mismatch on one barcode cannot map to another sample’s dual pair. Chimeric exact-piece pairs (BC1 of A + BC2 of B) stay unmatched under mm=1 in the 200k audit.

## How to re-run

```bash
seqmux demux -i sub_R1.fastq -I sub_R2.fastq \
  -b tests/fixtures/I464-469erdai_barcode_and_name.csv \
  -o /tmp/mm1 --counts-only --mismatches-1 1 --mismatches-2 1 -t 8 --force
python3 scripts/audit_mismatch1.py --table tests/fixtures/I464-469erdai_barcode_and_name.csv \
  --r1 sub_R1.fastq --r2 sub_R2.fastq --out audit.tsv
```
