# Full I464 Production Dataset mimalloc Validation Report

- **Dataset**: I464 Full Production (140,771,720 pairs)
- **Threads**: `-t 12`
- **Compression**: `--compression-level 1`
- **Build**: `target/release/seqmux` (`--features mimalloc-allocator`)
- **Correctness Gate**: **PASS**

## 1. Metrics Comparison

| Metric | Baseline v0.3.0 (glibc) | mimalloc-allocator | Delta / Speedup |
|:---|---:|---:|---:|
| **Wall Clock** | 406.80 s (6.78 min) | 195.23 s (3.25 min) | **2.08x** (+52.0%) |
| **Throughput** | 346,047 pairs/s | 721,056 pairs/s | **+375,009 pairs/s** |
| **Peak RSS** | 152.2 MB | 466.1 MB | +313.9 MB |
| **User CPU** | - | 1348.7 s | - |
| **Sys CPU** | - | 25.3 s | - |
| **CPU Utilization** | - | 703% | - |

## 2. Correctness Verification

| Metric | Expected (v0.3.0) | mimalloc Measured | Match |
|:---|---:|---:|:---:|
| Total Pairs | 140,771,720 | 140,771,720 | PASS |
| Assigned | 87,477,095 | 87,477,095 | PASS |
| Unassigned | 53,294,625 | 53,294,625 | PASS |
| Canonical | 45,393,500 | 45,393,500 | PASS |
| Swapped | 42,083,595 | 42,083,595 | PASS |

## 3. Engineering Conclusion
- 吞吐表现: 全量 1.41 亿 pairs 实际吞吐为 721,056 pairs/s，相比 baseline 提速 2.08x（耗时从 6分47秒 缩减至 3分15秒）。
- 内存开销: 峰值内存稳定在 466.1 MB，与 2M 采样测得的 467 MB 完全一致，证明内存严格有界（Bounded Memory），不会随 read 总数线性增长。
- 正确性: 全部统计指标与生产基准 100% 逐字吻合。
