# SeqMux 开发计划

> 工具名：`seqmux`  
> 目标平台：Linux x86_64 / Windows x86_64（MSVC）  
> 语言：Rust stable  
> 文档日期：2026-09-25

**当前状态：v0.3.0 RELEASED；v0.4 计划启动。** 生产配置与全量数字见 `docs/REAL_DATA.md`。第 2–8 节是到 v0.2.1 为止的历史记录，第 9 节是 v0.3.0 记录。v0.2.0、v0.2.1 与 v0.3.0 已打 tag 并正式发布。

早期 Ultraplex 重写草案（`ultraplex-rs`、Ultraplex CSV、`--three-prime-only`）已废弃，不再作为实现目标。

---

## 1. 产品定位

SeqMux 是单二进制、跨平台、低依赖的 FASTQ demultiplexer，面向实验室双 barcode / 单 barcode 扩增子文库：

1. 单端与 paired-end FASTQ（含 gzip）
2. SeqMux 样本表 CSV（header 必填）
3. Dual barcode：PE 两种插入方向；SE 为 R1 5′ + R1 3′
4. barcode 中 `N` 作为 UMI，写入 header `rbc:`
5. Hamming mismatch + tie → unassigned
6. quality / 3′ adapter trim、最短长度过滤
7. 多线程有序写出，无临时 FASTQ、无 pigz/SLURM/shell

原则不变：**正确性先于速度**。Demultiplex 最危险的失败是静默分错 sample。

---

## 2. v0.1 状态（已完成）

版本：`0.1.1` 功能基线 + `0.1.2` 真实数据方向修复。

| 能力 | 状态 |
|------|------|
| 单端 / 双端 FASTQ、gzip 读写 | 完成 |
| `seqmux validate` | 完成 |
| 样本表 CSV（含列名别名、忽略多余列） | 完成 |
| Dual / single barcode | 完成 |
| PE 两种插入方向（`--orientation both`，默认） | 0.1.2 |
| 方向规范化（output R1 = Barcode1 端） | 0.1.2 |
| UMI / mismatch / ambiguous | 完成 |
| quality + adapter trim、min length | 完成 |
| chunked multithread + ordered writer | 完成 |
| Linux + Windows CI、release 打包 workflow | 完成 |
| TSO `N/I` | 完成（实验室数据未用） |
| Ultraplex CSV 兼容 | **明确不做** |

v0.1 冻结范围仍有效：无 GUI、无 BAM/FASTA、无 UMI dedup、无 Python API、无分布式。

---

## 3. 真实实验数据实测（2026-08-13）

数据：`03_PCR_analysis/25_464469erdaimix`

| 项 | 值 |
|----|----|
| 文库 | GST-Sso7d 饱和突变，35 samples（A–E × 7 fragments） |
| 输入 | `E260702007_L01_464469erdaimix_{1,2}.fq.gz`，PE150 |
| 规模 | **140,771,720** pairs（约 15 GB gzip） |
| 实验室 Python demux | exact 8-mer，**两种方向都认**，match **62.14%** |
| 方向比例 | BC1@R1 45.4M（51.9% of matched）；BC2@R1 42.1M（48.1%） |
| Barcode1 ∩ Barcode2 | 空（方向不会把不同 sample 撞在一起） |

同实验室 `23_I395erdaimix`（60.5M pairs）同样约 47/53 两种方向，match 65.2%。这是 **扩增子 PE 的系统行为**，不是 I464 特例。

### 3.1 v0.1.1 在真实数据上的问题

对 200,000 pair 子样本（与全量 match rate 一致）：

| 策略 | 分配率 |
|------|--------|
| v0.1.1 仅 canonical（BC1@R1, BC2@R2） | **31.47%** |
| 实验室 Python（两种方向） | **62.43%** |
| 漏掉的可分配 reads | 约一半，全部是 swapped 方向 |

这是功能缺陷，不是性能问题：工具会“成功跑完”但丢掉 ~48% 本应入 sample 的 pair。

### 3.2 已在 0.1.2 修复

- 默认 `--orientation both`
- swapped hit 后 `mem::swap` R1/R2，使输出 R1 带 Barcode1（与样本表 F 引物端一致）
- `--orientation canonical` 可恢复旧行为，便于对照
- 真实 5-pair fixture：`tests/fixtures/i464_real_*.fastq`

SeqMux 规范化方向是 **Barcode1 → 输出 R1**。实验室 Python 脚本把 Barcode2 端写成输出 R1。下游若依赖旧脚本的 R1/R2 约定，用 `--no-canonicalize`。

### 3.3 v0.2 之前的观察（历史记录）

下列条目是 0.1.x 当时的笔记。除第 1 条（默认 mismatch=0，有意保留）和第 7 条（仍不把 round/fragment 写入 summary）外，都已在 v0.2 完成。

1. **~38% unmatched 仍在**（全量 53.3M pairs）。这是 exact match 的结果，与实验室 Python 一致；mismatch=1 的 QC 见 `docs/MISMATCH_QC.md`。
2. **已完成：** barcode 匹配在 quality/adapter trim 之前。
3. **已完成：** I464 / I395 全量已跑通，并与 Python 逐 sample 对照。v0.3 又在 gzip level 1 写出路径上复核，见第 9 节。
4. **已完成：** `--max-reads` / `--skip-reads`。
5. **已完成：** gzip 使用 `flate2` + `zlib-rs`。v0.2.0 counts-only 的墙钟见 `docs/REAL_DATA.md` 历史表；当前生产数字是 v0.3 gzip 写出。
6. **已完成：** stderr 进度含 pairs 与 pairs/s。
7. **仍不做：** 样本表里的 `library_round` / `PCR_product` / primers 不进入 summary。

---

## 4. v0.2 目标

主题：**实验室扩增子 demux 在真实数据上正确、可对照、可跑完全量。**

不做：GUI、UMI collapse、突变计数（仍由 `25_464469erdaimix` 下游脚本负责）、Ultraplex CSV 回归。

### 成功标准

- [x] I464 全量 1.41e8 pairs：assigned 与实验室 Python **逐 sample 计数一致**（exact、both orientation）
- [x] I395 全量同样对照通过
- [x] `--orientation canonical` 复现 v0.1.1 的 ~31% 分配率（对照用）
- [x] 1 mismatch 模式有明确 QC：多分配多少、是否把 chimeric pair 错分到其他 sample
- [x] barcode 匹配发生在 5′ quality trim 之前
- [x] 全量墙钟、peak RSS、reads/s 有记录（Linux）；内存不随 FASTQ 总大小增长
- [x] `--max-reads` 可用于冒烟
- [x] stderr 进度含 processed pairs 与 reads/s
- [x] README 写清方向约定 vs 实验室 Python 脚本差异
- [ ] Linux + Windows CI 仍全绿（本地 fmt/clippy/test 已通过；推送后看 Actions）

---

## 5. v0.2 Milestone

一次只做一个 milestone。每个结束后：

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all
cargo build --release
```

### HPC / 登录节点约束

开发机是 HPC 登录节点时，**重负载测试必须进计算队列，禁止在登录节点跑**。

| 环境 | 允许 |
|------|------|
| 登录节点 | 编辑代码、短单元测试、`cargo fmt` / `clippy`、小 fixture、`--max-reads` 冒烟 |
| 计算分区（如 `qcpu_23i`） | 真实数据 demux、全量 counts 对照、mismatch 审计、gzip 写盘 benchmark、长时间 / 高 CPU·IO 的 `cargo test`/`release` 以外重负载 |

提交方式：`sbatch` / `srun`。示例脚本：`tmp/run_m12_m14.slurm`；全量路径与对照约定见 `docs/REAL_DATA.md`。

### M10 — 真实数据回归脚手架 — **完成**

- `--max-reads N`、`--skip-reads`、`--counts-only`
- progress：pairs + k pairs/s + assigned %
- `docs/REAL_DATA.md`；summary `match_rate`

### M11 — 匹配顺序：barcode → trim — **完成**

顺序：match → canonicalize/trim barcode → quality → adapter → min-length。  
测试：`barcode_match_survives_low_quality_five_prime`。

### M12 — I464 / I395 全量对照 — **完成**

SLURM job 2313038（`qcpu_23i`，16 threads）。与实验室 Python exact demux **逐 sample / orientation / unassigned 完全一致**。差异见 `docs/COMPATIBILITY.md`。

### M13 — mismatch=1 安全性评估 — **完成**

默认保持 0。I464 全量 mm=1 多分配 +3.33M（+2.36 pp），无 ambiguous、无 sample 计数下降；子样本 chimeric 0 吸收。见 `docs/MISMATCH_QC.md`。

### M14 — I/O 与进度 — **完成**（v0.2.0 历史测量）

v0.2.0 counts-only 全量约 390 k pairs/s，RSS 约 16 MB。该数字只描述当时的 counts-only 运行，不是 v0.3 生产吞吐。v0.3 在 gzip level 1 写出下的全量数字见第 9 节。已启用 `flate2` `zlib-rs` backend、thin LTO。未上 SIMD/unsafe。

### M15 — v0.2 发布 — **完成**

- tag `v0.2.0` 与 `v0.2.1` 已存在
- 可选 musl / aarch64，非必须

---

## 6. 建议 CLI 增量（v0.2）

已有（0.1.2）：

```text
--orientation both|canonical|swapped
--no-canonicalize
--mismatches-1 / --mismatches-2
```

v0.2 新增：

```text
--max-reads <INT>
--skip-reads <INT>
--counts-only
```

不要为了和实验室 Python 脚本参数名一致而加入 `pigz-threads` 等无意义项。

---

## 7. 架构约束（继续有效）

```text
Rust stable，单 Cargo package
byte-oriented FASTQ
compiled barcode + Hamming（禁止 5^L 字典）
chunked bounded pipeline，N workers，1 ordered writer
直接写最终文件，无 _tmp_thread_* ，无 cat/pigz/mv
Path/PathBuf，Linux + Windows CI
默认不覆盖已有输出（--force）
```

Barcode 匹配：mismatch=0 可用 HashMap；有 mismatch 时 Hamming。barcode 很短，先正确再 packed Hamming。

---

## 8. 测试策略（v0.2）

```text
Unit            matcher 方向 / canonicalize / trim 顺序
Integration     tests/fixtures/i464_real_*.fastq
Real-data       I464 / I395 全量 vs Python counts（仓库外数据）
Mismatch study  错分审计表
```

全量 FASTQ **不进 git**。文档只记录绝对路径与期望计数。  
Real-data / mismatch / benchmark 等重负载遵守上文「HPC / 登录节点约束」。

---

## 9. v0.3.0 RELEASED

功能与全量验证已完成，v0.3.0 已打 tag 并正式发布。

- 双线程并发解压 R1/R2（`-t > 1`）；`-t 1` 仍走同步 reader
- stride-1 adapter：固定 8-mer 查找表 + 短 overlap 4-mer 候选；generic DP 仍是参照
- 生产配置：`-t 12`，`--compression-level 1`（CLI 默认压缩级别仍是 6），mismatch 0，adapter 默认开启，orientation `both` 且 canonicalize
- I464 全量 346,047 pairs/s（406.80 s）；I395 全量 341,880 pairs/s（177.03 s）
- I464 35 个样本、I395 18 个样本与实验室 Python exact demux 一致

仍不做：GUI、UMI collapse、writer 分片、额外 reader、SIMD。

### Ultraplex benchmark（已完成对照，见 `docs/ULTRAPLEX_BENCH.md`）

I464 200k pairs、exact / `-q 0`：SeqMux 与 Ultraplex **逐 sample 计数一致**（canonical 与 both）。  
墙钟约 **10×（8 线程）～36×（1 线程）** 快于 Ultraplex；Ultraplex 必须 `--dont_build_reference`，且勿用本机版 `--ignore_no_match`（TypeError）。

---

# SeqMux v0.4 Development Plan

## 0. 当前基线

当前正式版本：

```text
SeqMux v0.3.0
```

v0.3.0 已完成：

- single-end / paired-end FASTQ demultiplex
- single / dual barcode
- PE canonical + swapped orientation
- swapped mate canonicalization
- exact 8-bp dual-barcode fast path
- mismatch matching
- barcode `N` / UMI
- quality trimming
- 3′ adapter trimming
- concurrent PE decompression
- bounded multi-worker processing pipeline
- ordered writer
- gzip output
- counts-only
- FASTQ validate
- Linux + Windows CI / release binary
- I464 / I395 全量真实数据验证
- 与实验室 Python exact demux 逐 sample 计数一致

v0.3.0 当前生产参考：

```text
I464
140,771,720 pairs
~346k pairs/s
35 samples exact parity

I395
60,523,000 pairs
~342k pairs/s
18 samples exact parity
```

核心原则继续冻结：

> Correctness > performance。Demultiplex 最危险的问题是 silent mis-assignment，而不是少跑 5–10% 的速度。

---

## 1. v0.4 总目标

v0.4 不进行大规模功能扩张。

本阶段只解决以下四件事：

```text
1. 修复发布状态/documentation drift
2. 增强 mismatch > 0 时的 assignment ambiguity protection
3. 建立 SeqMux vs fqtk 公平 benchmark
4. 建立长期可重复的 performance regression framework
5. 为 crates.io / Bioconda 分发做好准备
```

不把以下内容纳入 v0.4 主线：

```text
❌ fqtk-style arbitrary read structures
❌ I1/I2 任意 index-read protocol
❌ cellular barcode engine
❌ UMI collapse/dedup
❌ BAM support
❌ GUI
❌ SIMD
❌ unsafe optimization
❌ writer sharding
❌ distributed processing
❌ Python API
```

除非 benchmark 明确证明当前 architecture 存在生产瓶颈，否则不要为了“架构更先进”重写 pipeline。

---

## M16 — v0.3.0 发布状态收口

### 目标

修复仓库中已经与 GitHub release 状态不一致的文档。

当前事实：

```text
Cargo.toml = 0.3.0
GitHub release v0.3.0 已存在
```

但部分文档仍写：

```text
Release candidate
Not tagged
等待 tag / push
```

这些必须统一。

### 修改范围

检查至少：

```text
CHANGELOG.md
DEVELOPMENT_PLAN.md
README.md
docs/REAL_DATA.md
```

#### CHANGELOG

将：

```text
0.3.0 — Release candidate. Not tagged.
```

改为正式 release 状态。

不要修改已经冻结的 benchmark 数字。

#### DEVELOPMENT_PLAN

将：

```text
v0.3.0 release candidate
等待用户审阅
不要 tag/push
```

改为：

```text
v0.3.0 RELEASED
```

然后追加本 v0.4 计划。

### Acceptance

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all
```

必须通过。

文档修改不得改变程序行为。

---

## M17 — Matching ambiguity hardening

**状态：COMPLETE；local gates PASS（110 tests）；production regression PASS。**
证据归档在 `benchmarks/m17_regression/`（文本；FASTQ 输出不入 git）。
执行过程中发现 `--min-mismatch-delta` 在现有 I464/I395 barcode 表上是惰性参数，
原因与含义见 `benchmarks/m17_regression/README.md`。

### 目的

SeqMux 当前：

```text
best distance wins
exact tie → ambiguous/unassigned
```

对于：

```text
mismatch = 0
```

该行为已经足够安全。

但当用户启用：

```text
--mismatches-1 1
--mismatches-2 1
```

时，目前只防止完全 tie，没有限制：

```text
best = 1
second-best = 2
```

这种低 margin assignment。

参考 fqtk，引入：

```text
--min-mismatch-delta
```

但必须保持 v0.3 backward compatibility。

---

### 语义定义

定义：

```text
best_score = 最优 eligible sample 的总 mismatch
second_best_score = 次优不同 eligible sample 的总 mismatch

delta = second_best_score - best_score
```

candidate 首先必须分别满足：

```text
d1 <= mismatches_1

以及 dual 模式：

d2 <= mismatches_2
```

然后：

```text
score = d1 + d2
```

single barcode：

```text
score = d1
```

Dual PE 先取每个 sample 在允许 orientation 下的最小 score，再比较不同
samples。同 sample 的 canonical/swapped 同分时保留 canonical（仅允许
swapped 时保留 swapped），不据此判 Ambiguous。没有第二个 eligible sample
时直接 assignment；threshold 外的 sample 不参与 margin。

assignment 条件：

```text
best candidate exists
AND
delta >= min_mismatch_delta
```

---

### Backward compatibility

CLI 新增：

```text
--min-mismatch-delta <INT>
```

默认（兼容 v0.3.0）：

```text
0
```

但必须保留当前 tie 行为。

因此实现时注意：

```text
best == second-best
```

无论 delta 参数如何，都必须继续：

```text
Ambiguous
```

也就是说不能因为默认 delta=0 而让 tie 被 assignment。

建议逻辑：

```text
if no best:
    NoMatch

if multiple different samples share best score:
    Ambiguous

if second best from a different eligible sample exists
   AND second_best - best < min_mismatch_delta:
    Ambiguous

otherwise:
    Assigned
```

这样：

```text
--min-mismatch-delta 0
```

与 v0.3 行为完全一致。

---

### Exact fast path

以下生产路径不得出现性能回退：

```text
8 bp + 8 bp
ACGT only
mismatch1 = 0
mismatch2 = 0
```

应继续使用：

```text
ExactDual8Matcher
```

不要为了支持 delta 强制 exact path 回退到 generic matcher。

因为 exact lookup 的 key 唯一性在 config parsing 阶段已经验证。

---

### MatchResult

如有必要，可以增强内部结构，例如：

```rust
Assigned {
    sample_idx,
    orientation,
    mismatches,
}
```

但不要在没有必要的情况下修改公开输出格式。

summary 可以考虑增加：

```text
ambiguous
```

如果已经有该字段，则保持格式兼容。

不要随意增加大量 debug columns。

---

### 测试

必须覆盖：

#### single barcode

```text
best=0, second=2, delta=2 → assign
best=1, second=2, min_delta=2 → ambiguous
best=1, second=3, min_delta=2 → assign
tie 1/1 → ambiguous
```

#### dual barcode

同时验证：

```text
d1+d2 score
per-barcode mismatch threshold
combined second-best delta
```

例如：

```text
sample A: d1=1 d2=0 total=1
sample B: d1=1 d2=1 total=2
min_delta=2

→ ambiguous
```

#### orientation

必须覆盖：

```text
canonical candidate
swapped candidate
canonical vs swapped competing candidate
```

尤其防止：

```text
different samples
canonical / swapped equal score
```

被错误静默 assignment。

同 sample canonical/swapped equal score 必须继续 Match，保持既有 canonical
优先行为；不能把这两个 orientation 当作 best/second-best samples。

#### N / UMI

不能破坏已有 `N` semantics。

#### property tests

扩展现有 `proptest`。

至少检查：

```text
result independent of sample-table row order
exact-match unique barcode always stable
increasing min-mismatch-delta cannot create additional assignments
```

最后一条非常重要：

```text
Assigned(delta=N+1)
⊆
Assigned(delta=N)
```

---

### Real-data regression

重新跑：

```text
I464 mismatch=0
I395 mismatch=0
```

必须完全保持：

```text
total
assigned
unassigned
canonical
swapped
per-sample counts
```

与 v0.3 production baseline 一致。

然后重新检查：

```text
I464 mismatch=1
```

至少比较：

```text
delta=0
delta=1
delta=2
```

记录：

```text
assigned
ambiguous
unassigned
per-sample changes
assignment gained/lost
```

目标不是证明 delta=2 一定最好。

目标是得到数据。

---

### M17 完成标准

必须满足：

```text
v0.3 exact behavior unchanged
exact 8+8 fast path retained
new ambiguity semantics deterministic
row-order invariant
tests green
I464/I395 mismatch=0 exact regression PASS
```

不得在本 milestone 做性能优化。

---

## M18 — SeqMux vs fqtk benchmark

**状态：MEASURED / CORRECTNESS GATE PASS；报告已写，矩阵有缺口未补。**
结果与完整记录见 `docs/FQTK_BENCHMARK.md`；原始机器可读数据
`benchmarks/fqtk/benchmark_results.tsv`（36 行 measured + warm-up）。

### 已确立

- 2M 与全量 140,771,720 pairs：SeqMux vs fqtk **逐样本计数一致**，
  且 2M 做过流式 FASTQ 内容比对（read ID、顺序、去 barcode 后序列）。
  见 `benchmarks/fqtk/provenance/gate_result.txt`、`full_parity_result.txt`。
- t12/c1：fqtk 比 SeqMux 快 **1.37×**（2M median 4.22 s vs 5.77 s）；
  全量 1.34×（299.25 s vs 401.36 s，单次运行）。
- SeqMux peak RSS 低于 fqtk（全量 138 MB vs 1331 MB）。
- c6 下 SeqMux 慢约 3.8×（16.0 s vs 4.2 s）。

### 两个必须随数字一起引用的 caveat

1. **nominal compression level 不等于等量压缩工作。** 同为 level 1，SeqMux 写出
   225.9 MB，fqtk 150.4 MB（+50%）；而 level 6 下 SeqMux 反而更小
   （88.6 MB vs 97.9 MB）。因此 level-1 的墙钟对比不是"两个工具做同样的压缩"。
   生产链路取 `--compression-level 1`（I464 产出 10.9 GB），此处有真实存储/IO 成本。
   根因分析**不在 M18**，留给 M21，且需先 profile 证明 writer 是瓶颈。
2. **8–16 线程内无扩展性**（5.81/5.77/5.92 s）。该 workload 只有 5 秒，处于
   startup/IO 噪声区。

### M18 补齐与收口状态（已完成）

- 线程矩阵缺口已补齐：SeqMux 在 t1 与 t4（c1、c6）完成 3 次重复扩展测量；fqtk 因强制 `-t >= 5`，在 t5 与 t6 提供了最低并发对齐基准。
- I395 2M 基准已完成：correctness gate 100% 逐样本及 FASTQ 内容对齐（18 样本，610,058 assigned）；完成 t1/t4/t5/t8/t12/t16 矩阵测试（3 重复）。
- 全量运行保持单次（作为生产级全量对照，豁免 3 重复）。
- M18 已达标并标定 **PASS**。

### 原则

这是 scientific benchmark，不是 marketing benchmark。

不要预设 SeqMux 应该更快。

最终允许的结论只有数据支持的内容。

---

### 2.1 固定版本

记录：

```text
SeqMux version
SeqMux commit SHA

fqtk version
fqtk commit SHA

rustc --version
CPU model
CPU count
RAM
OS
filesystem
```

fqtk 使用固定 release/tag。

不要 benchmark 一个不断变化的 `main`。

---

### 2.2 Primary benchmark：canonical exact demux

这是最公平的 apples-to-apples benchmark。

使用：

```text
I464
```

先用：

```text
2M pairs fixed subset
```

完成 benchmark matrix。

然后在必要时运行：

```text
140.8M full dataset
```

---

#### SeqMux semantic configuration

类似：

```bash
seqmux demux \
  -i R1.fq.gz \
  -I R2.fq.gz \
  -b barcodes.csv \
  -o OUT \
  --orientation canonical \
  --no-adapter \
  -q 0 \
  --mismatches-1 0 \
  --mismatches-2 0 \
  --compression-level X \
  -t N
```

确认无 quality trimming。

Barcode 默认 trim。

---

#### fqtk semantic mapping

构造 metadata：

```text
sample_id
barcode
```

其中：

```text
barcode = Barcode1 + Barcode2
```

read structure：

```text
R1: 8B+T
R2: 8B+T
```

例如：

```bash
fqtk demux \
  --inputs R1.fq.gz R2.fq.gz \
  --read-structures 8B+T 8B+T \
  --sample-metadata metadata.tsv \
  --output OUT \
  --max-mismatches 0 \
  --header-format unmodified \
  --threads N \
  --compression-level X
```

必须首先用小 fixture 验证 semantics。

---

### 2.3 Correctness gate BEFORE performance

在记录任何性能数据以前，必须证明：

```text
SeqMux assigned count == fqtk assigned count
SeqMux unassigned == fqtk unmatched
per-sample counts identical
```

至少抽查输出：

```text
record counts
read IDs
template sequence after barcode removal
paired synchronization
```

如果 semantics 不一致：

> 停止 benchmark，不得比较速度。

---

### 2.4 Benchmark matrix

第一轮：

```text
threads:
1
4
8
12
16
```

compression：

```text
1
6
```

每种组合：

```text
3 repetitions
```

报告：

```text
median wall
min/max wall
pairs/s
user CPU
sys CPU
CPU %
peak RSS
output bytes
```

不要只报最快一次。

---

### 2.5 Warm-up / cache

每个工具至少先：

```text
1 warm-up run
```

正式重复不要把 warm-up 计入结果。

如果 filesystem cache 不能严格控制，要在报告中明确。

不要假装这是 cold-cache benchmark。

---

### 2.6 Secondary benchmark：SeqMux native both-orientation

SeqMux：

```text
--orientation both
```

这是 SeqMux 的核心实验室场景。

fqtk 没有完全等价的“一 sample 两种 orientation”原生语义。

允许构造：

```text
sample_A_canonical
sample_A_swapped
```

两组 metadata rows，然后结果 aggregate 回：

```text
sample_A
```

但必须明确标为：

```text
protocol emulation benchmark
```

不能和 primary benchmark 混为一个 apples-to-apples benchmark。

---

### 2.7 不公平 benchmark 禁止项

禁止：

```text
SeqMux --counts-only
vs
fqtk full FASTQ output
```

禁止：

```text
SeqMux --no-adapter
vs
另一个工具做 adapter trimming
```

禁止：

```text
不同 compression level
```

然后宣称工具速度差异。

---

### 2.8 Benchmark artifact

新增例如：

```text
docs/FQTK_BENCHMARK.md
scripts/run_fqtk_benchmark.py
scripts/compare_fqtk_counts.py
```

建议 benchmark script 产生机器可读：

```text
benchmark.tsv
```

字段至少：

```text
tool
version
commit
dataset
pairs
mode
threads
compression_level
replicate
wall_seconds
pairs_per_second
user_seconds
sys_seconds
cpu_percent
peak_rss_kb
output_bytes
assigned
unassigned
```

最终文档从该 TSV 汇总。

不要手工抄 benchmark 数字。

---

## M19 — Performance regression framework

目前 SeqMux 已有真实 benchmark，但性能记录主要存在：

```text
docs/
scripts/
历史 SLURM output
```

v0.4 要把它变成可重复工程设施。

---

### 目标

新增：

```text
scripts/run_perf_regression.py
```

或者扩展现有 production validation script。

它应：

1. 接收 FASTQ / barcode table
2. 固定 command configuration
3. 运行 release binary
4. 捕获 `/usr/bin/time -v` 或等价信息
5. 解析 summary
6. 输出 machine-readable TSV/JSON
7. 与 optional baseline 比较

---

### 两级 benchmark

#### Level A — repository smoke

使用小型 deterministic fixture。

目标：

```text
correctness
pipeline invariants
command reproducibility
```

不在 GitHub Actions 中强制 wall-clock threshold。

共享 CI 性能噪声太大。

---

#### Level B — HPC production benchmark

使用：

```text
I464 2M
I464 full
I395 full
```

真实 FASTQ 不进入 git。

配置路径通过 CLI / environment 提供。

---

### Baseline manifest

建议加入类似：

```text
benchmarks/baselines/v0.3.0.tsv
```

只记录已经冻结的生产结果和配置，不记录私有 FASTQ 内容。

至少包含：

```text
dataset identifier
record count
SeqMux version
threads
compression
adapter mode
mismatch
orientation
wall
RSS
assigned
unassigned
```

不要把机器相关 wall time 当作跨机器 universal threshold。

---

### Regression policy

正确性：

```text
任何 count 变化 → FAIL
```

性能：

```text
>10% median slowdown
```

先报告：

```text
WARNING / INVESTIGATE
```

不要直接在公共 CI fail。

只有固定 benchmark host 上才能使用硬 gate。

---

## M20 — Distribution readiness

目标不是立即发布，而是让 SeqMux 达到可发布条件。

---

### 4.1 crates.io readiness

检查：

```text
Cargo.toml metadata
README links
license
repository
keywords
categories
MSRV
```

运行：

```bash
cargo package
cargo package --list
```

确认：

```text
不包含真实 FASTQ
不包含 tmp/
不包含 benchmark outputs
不包含 SLURM logs
不包含私有路径文件
```

然后测试：

```bash
cargo install --path .
```

以及从 package artifact 本地安装。

实际：

```text
cargo publish
```

必须等待用户明确授权。

---

### 4.2 Release binary

继续保持：

```text
Linux x86_64
Windows x86_64
```

检查 release workflow：

```text
version
tag
binary --version
archive naming
SHA256
```

建议 release asset 增加 checksums 文件：

```text
SHA256SUMS
```

---

### 4.3 Bioconda preparation

准备：

```text
meta.yaml
build.sh
```

或者独立 bioconda-recipe PR 所需内容。

只准备 recipe 和本地验证。

不要未经用户许可直接提交外部 PR。

---

## M21 — Conditional performance work

这是条件 milestone。

M18 benchmark 完成前：

> 禁止开始 M21。

只有 profile / benchmark 明确证明瓶颈后才允许实施。

---

### 5.1 如果 writer/compression 是瓶颈

首先 profile：

```text
writer CPU
gzip CPU
worker idle fraction
channel blocking
```

优先尝试低复杂度优化。

例如评估：

```text
compression worker pool
```

但必须保持：

```text
per-sample read order
pair synchronization
deterministic assignment
bounded memory
output failure safety
```

不要直接复制 fqtk pooled-writer architecture。

先验证 SeqMux 是否真的需要。

---

### 5.2 如果 barcode matching 是瓶颈

优先：

```text
retain exact packed fast path
optimize generic matcher only with profile evidence
```

禁止：

```text
SIMD for SIMD's sake
unsafe for benchmark numbers
```

---

### 5.3 如果 adapter trimming 是瓶颈

已有：

```text
packed seed fast path
generic DP correctness reference
differential tests
adversarial tests
```

任何优化必须继续 differential-check against reference implementation。

---

## M22 — v0.4 release gate

只有前面 milestone 完成后再进入。

最低要求：

```text
M16 PASS
M17 PASS
M18 PASS
M19 PASS
M20 readiness PASS
```

M21：

```text
optional
```

不是 release blocker。

---

### Release correctness gates

必须：

```bash
cargo fmt --check

cargo clippy \
  --all-targets \
  --all-features \
  -- -D warnings

cargo test --all

cargo build --release
```

Linux + Windows CI 全绿。

---

### Production regression

至少重新跑：

```text
I464 mismatch=0
I395 mismatch=0
```

必须：

```text
total exact
assigned exact
unassigned exact
per-sample exact
orientation canonical exact
orientation swapped exact
```

---

### Output safety

再次验证：

```text
default refuses overwrite
--force works
truncated / invalid FASTQ fails closed
R1/R2 mismatch fails
worker panic propagates
missing chunk invariant fails
gzip finalize succeeds
```

---

## 7. Architecture invariants

以下约束作为 v0.4 architecture contract：

```text
Reader
   ↓
bounded chunks
   ↓
N workers
   ↓
OrderedWriter
   ↓
per-sample FASTQ
```

必须保持：

```text
bounded memory
input-order preservation
paired synchronization
single authoritative assignment
workers do not own persistent sample output files
no temp FASTQ merge workflow
```

---

## 8. Scientific correctness invariants

必须保持：

### barcode

```text
match original sequence before quality trimming
```

### paired end

```text
R1/R2 ID must remain synchronized
```

### orientation

默认：

```text
canonical + swapped
```

canonicalization：

```text
output R1 = Barcode1 side
```

### ambiguity

```text
never silently select equal best candidates from different samples
```

### exact mode

```text
mismatch=0 must reproduce v0.3 production counts
```

---

## 9. 明确不做的 feature creep

Gemini 不得自行加入：

```text
fqtk read-structure DSL
I1/I2 generalized indexing
10x barcode support
BAM
UMI consensus
deduplication
distributed mode
GPU
SIMD
unsafe
GUI
web server
```

如果发现这些可能有价值，只记录：

```text
Future consideration
```

不要实现。

---

## 10. 工作方式

一次只执行一个 milestone。

顺序严格：

```text
M16
 ↓
M17
 ↓
M18
 ↓
M19
 ↓
M20
 ↓
M21 only if evidence justifies it
 ↓
M22
```

每个 milestone 完成后必须汇报：

```text
1. 修改了什么
2. 为什么修改
3. 哪些文件改变
4. tests
5. correctness evidence
6. performance evidence（如适用）
7. unresolved risks
8. READY / NOT READY for next milestone
```

---

## 11. Git 操作约束

Gemini 可以：

```text
修改代码
运行测试
运行 benchmark
生成报告
```

未经用户明确要求不要：

```text
push
创建 GitHub release
cargo publish
提交 Bioconda PR
```

commit 是否创建由用户当前指令决定。

不要自动改变已有 release tag。

---

## 12. 当前状态与下一条任务

```text
M16 — v0.3.0 发布状态收口已完成。
M17A — COMPLETE & CLOSED (Commit 2d9e5ba)。严谨 A/B 评测与全量 1.41 亿 pairs 验证完成：拒绝 fat LTO，保留 mimalloc 作为可选 feature。
M17B — 重构为真正置信度余量语义（全局近邻不同样本竞争，定向合并，全量测试 107 tests 全部通过，真实数据几何证明与扫描完成）。
M18 — PASS；correctness gate PASS；报告 docs/FQTK_BENCHMARK.md；矩阵缺口已补齐（t1/t4 与 I395 矩阵）。
M19 — PASS；两级回归框架已落地（scripts/run_perf_regression.py 与 benchmarks/baselines/v0.3.0.tsv）；
      Level A (smoke) 与 Level B (i464_2m, i464_full, i395_full) 全量端到端验证通过，报告归档 benchmarks/perf_regression_report.tsv。
M20 — 就绪，待启动（Distribution readiness：crates.io readiness, release binary, Bioconda）。
```

### 仓库卫生决策（M18 收口时确定）

- `.gitignore` 新增 `benchmarks/**/runs*/`、`benchmarks/**/gate_*`、`*.fq.gz`、
  `*.fastq.gz`。这些目录单是 I464 全量就有 21 GB，绝不能进索引。
- `benchmarks/` 下**保留**可复现的小记录：`benchmark_results.tsv`、
  `provenance/`、`commands.log`、`i464_metadata.tsv`、`slurm_*.out/err`。
- `benchmarks/m17_regression/**` 用 `!` 例外强行跟踪（含 `*.summary.tsv`），
  因为它是 M17 唯一的生产回归证据；`tmp/` 被忽略后这些证据会丢失。
- 提交前用 `git add -An` 确认零 FASTQ 进入索引。

---

## 参考

- 本仓库 README / CHANGELOG
- 实验室对照：`/hpcfs/fhome/caizhh/03_PCR_analysis/25_464469erdaimix`
- 实验室对照：`/hpcfs/fhome/caizhh/03_PCR_analysis/23_I395erdaimix`
- 旧 Ultraplex 仅作历史灵感：https://github.com/ulelab/ultraplex
