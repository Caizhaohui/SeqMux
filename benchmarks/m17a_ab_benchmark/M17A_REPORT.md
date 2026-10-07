# SeqMux M17A Benchmark Report — Release Build / Allocator Strict A/B

## 1. Environment & Hardware Frozen Provenance

- **Host**: `bnode1.tibhpc.net` (SLURM jobs `2893788` and `2893854`)
- **CPUs allocated**: `16` (dedicated allocation, `Cpus_allowed_list` verified)
- **Rust compiler**: `rustc 1.83.0 (90b35a623 2024-11-26)` / `cargo 1.83.0`
- **Methodology**: 1 warm-up + 5 interleaved measured repetitions across 4 variants (Round 1: A B C D, Round 2: D C B A, Round 3: B D A C, Round 4: C A D B, Round 5: A C B D).
- **Provenance files**: Archived in `benchmarks/m17a_ab_benchmark/provenance/` (git commit, rustc/cargo versions, lscpu, slurm env).

## 2. 2M Production Workload Results (Median of 5 Replicates)

Dataset: I464 2,000,000 paired-read benchmark dataset (`tmp/upx_bench/i464_2000000_*.fq.gz`)
Command: `-t 12 --compression-level 1`

| Variant | Description | Wall (s) | Min–Max (s) | pairs/s | Speedup | User CPU | Sys CPU | Peak RSS (MB) |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|
| **A** | A_baseline (thin LTO, codegen-units=16, default glibc) | 6.01 | 5.86–6.09 | 332,779 | **1.00x** | 35.60s | 16.43s | 125.7 |
| **B** | B_fat_lto (fat LTO, codegen-units=1, default glibc) | 6.14 | 5.86–6.35 | 325,733 | **0.98x** | 36.89s | 18.23s | 146.1 |
| **C** | C_mimalloc (thin LTO, codegen-units=16, mimalloc) | 3.54 | 3.48–4.54 | 564,972 | **1.70x** | 18.51s | 0.55s | 467.1 |
| **D** | D_both (fat LTO, codegen-units=1, mimalloc) | 3.60 | 3.40–3.79 | 555,556 | **1.67x** | 18.43s | 0.55s | 457.6 |

## 3. Diagnostic Counts-Only Workload & Build Trade-Offs

Command: `-t 12 --counts-only`

| Variant | Diagnostic Wall (s) | Diag Speedup | Clean Build Time (s) | Binary Size (MB) |
|:---|---:|---:|---:|---:|
| **A** | 5.61 | 1.00x | 10.5s | 3.62 MB |
| **B** | 5.65 | 0.99x | 36.0s | 3.00 MB |
| **C** | 1.96 | 2.86x | 7.5s | 3.76 MB |
| **D** | 2.02 | 2.78x | 29.0s | 3.15 MB |

## 4. Full Production I464 Validation (140,771,720 pairs)

Dataset: Full I464 (`E260702007_L01_464469erdaimix_*.fq.gz`, 15.1 GB compressed FASTQ)
Command: `-t 12 --compression-level 1` on 16 CPUs (`qcpu_18i`, Job `2893854`)

| Metric | Baseline v0.3.0 (glibc) | mimalloc-allocator | Delta / Speedup |
|:---|---:|---:|---:|
| **Wall Clock** | 406.80 s (6.78 min) | 195.23 s (3.25 min) | **2.08x** (+52.0% time reduction) |
| **Throughput** | 346,047 pairs/s | **721,056 pairs/s** | **+375,009 pairs/s** |
| **Peak RSS** | 152.2 MB | **466.1 MB** | +313.9 MB (bounded, identical to 2M) |
| **User CPU** | - | 1,348.7 s | - |
| **Sys CPU** | - | 25.3 s | (lock contention negligible) |

## 5. Correctness Gate Verification

- **2M Dataset Parity**: All 4 variants produced identical sample counts (1,251,876 assigned, 748,124 unassigned, 629,825 canonical, 622,051 swapped, 16,542 trimmed). Decompressed FASTQ SHA-256 for all 70 sample streams matched Baseline A byte-for-byte.
- **Full 140M Dataset Parity**: 140,771,720 pairs processed. Exactly 87,477,095 assigned, 53,294,625 unassigned, 45,393,500 canonical, 42,083,595 swapped. 100% matched production v0.3.0 baseline across all 35 samples.
- **Test Suites**: `cargo test --all` and `cargo test --all --features mimalloc-allocator` both passed 100% (105 integration/unit tests). `cargo clippy --all-targets --all-features` produced 0 warnings.

## 6. M17A Final Verdict

```text
fat LTO:              REJECT
                      - Performance: 0.98x (no improvement / marginal slowdown)
                      - Compilation cost: 36.0s vs 10.5s (3.4x longer)
                      - Decision: Do not adopt fat LTO or release-fat profile.

mimalloc:             KEEP AS OPTIONAL CARGO FEATURE (`--features mimalloc-allocator`)
                      - Performance: 1.70x on 2M, 2.08x on full 140M (721k pairs/s), 2.86x on counts-only (1.02M pairs/s)
                      - Lock contention: Sys CPU collapsed from ~16s down to 0.5s on 2M, and 25s on full 140M
                      - Memory trade-off: RSS increases from 152 MB to 466 MB, but is strictly bounded across 140M pairs
                      - Decision: Keep as an optional feature for HPC high-throughput deployments.
                                  Default SeqMux release builds continue using the system allocator for minimal RSS footprint.

RECOMMENDED PRODUCTION BUILDS:
- Default Release: cargo build --release (thin LTO, codegen-units=16, default allocator, RSS ~150 MB)
- HPC High-Throughput: cargo build --release --features mimalloc-allocator (throughput ~720k pairs/s, RSS ~466 MB)
```
