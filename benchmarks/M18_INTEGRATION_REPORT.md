# SeqMux M18 — v0.4.0 Integration & Regression Report

- **Date**: 2026-10-07 16:38:44
- **Host**: `bnode2`
- **Test Dataset**: I464 2,000,000 pairs & I395 2,000,000 pairs
- **Workload**: `-t 12 --compression-level 1` (production gzip) / `--counts-only` (matrix)

## 1. M18.2 Production Behavior Regression (I464 2M, mismatch=0, delta=0)

| Build Variant | Wall (s) | Throughput (pairs/s) | Speedup | User CPU | Sys CPU | CPU % | Peak RSS (MB) | Output Bytes |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| **Standard Allocator** | 5.77s | 346,620 | 1.00x | 35.26s | 14.68s | 842% | 115.9 MB | 194,708,499 B |
| **mimalloc Allocator** | 3.50s | 571,429 | **1.65x** | 18.66s | 0.52s | 548% | 452.5 MB | 194,708,499 B |

### Output Verification against v0.3.0 Baseline:
- Total reads: 2,000,000
- Assigned: 1,251,876 (62.59%)
- Unassigned (includes ambiguous): 748,124 (37.41%)
- Ambiguous subset: 0 (0.00%)
- Orientation Canonical: 629,825
- Orientation Swapped: 622,051
- Adapter Trimmed: 16,542

## 2. M18.3 Output Parity Gate (Decompressed FASTQ SHA-256)
- **Files checked**: 72 sample files (R1 & R2)
- **Result**: **100% BYTE-FOR-BYTE IDENTICAL** across Standard and mimalloc builds.

## 3. M18.4 Mismatch-Delta Behavior Matrix

### I464 2M (mismatch=1):
| Delta | Allocator | Assigned | Unassigned (total) | Ambiguous subset | Status |
|:---|:---|---:|---:|---:|:---|
| delta=0 | Standard | 1,308,507 | 691,493 | 0 | Baseline (no delta) |
| delta=0 | mimalloc | 1,308,507 | 691,493 | 0 | Identical |
| delta=1 | Standard | 1,308,507 | 691,493 | 0 | Geometry lower bound delta>=4 satisfied |
| delta=2 | Standard | 1,308,507 | 691,493 | 0 | Geometry lower bound delta>=4 satisfied |
| delta=5 | Standard | 1,308,421 | 691,579 | 86 | Confidence margin active (86 ambiguous) |
| delta=5 | mimalloc | 1,308,421 | 691,579 | 86 | Identical to Standard |

### I395 2M (mismatch=1):
| Delta | Allocator | Assigned | Unassigned (total) | Ambiguous subset | Status |
|:---|:---|---:|---:|---:|:---|
| delta=0 | Standard | 1,340,210 | 659,790 | 0 | Baseline |
| delta=2 | Standard | 1,340,210 | 659,790 | 0 | Geometry lower bound delta>=3 satisfied |
| delta=4 | Standard | 1,340,209 | 659,791 | 1 | Confidence margin active (1 ambiguous) |
| delta=4 | mimalloc | 1,340,209 | 659,791 | 1 | Identical to Standard |
