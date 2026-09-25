# SeqMux 开发计划

> 工具名：`seqmux`  
> 目标平台：Linux x86_64 / Windows x86_64（MSVC）  
> 语言：Rust stable  
> 文档日期：2026-09-25

**当前状态：v0.3.0 release candidate。** 生产配置与全量数字见 `docs/REAL_DATA.md`。第 2–8 节是到 v0.2.1 为止的历史记录，不是未完成工作。v0.2.0 与 v0.2.1 已打 tag。

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

## 9. v0.3.0 release candidate

功能与全量验证已完成。在用户明确要求之前不打 tag、不 push。

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

## 10. 给 coding agent 的下一条任务

```text
v0.3.0 生产验证已完成，等待用户审阅报告。

不要在用户明确要求之前 commit、tag 或 push。
不要再做微优化，除非用户指出全量 I464/I395 运行中的新瓶颈。
```

---

## 参考

- 本仓库 README / CHANGELOG
- 实验室对照：`/hpcfs/fhome/caizhh/03_PCR_analysis/25_464469erdaimix`
- 实验室对照：`/hpcfs/fhome/caizhh/03_PCR_analysis/23_I395erdaimix`
- 旧 Ultraplex 仅作历史灵感：https://github.com/ulelab/ultraplex
