# SeqMux 开发计划 (from Ultraplex Rust rewrite plan)

> 工作名：`ultraplex-rs`（仅作为开发阶段占位名，后续可更名）  
> 目标平台：Linux x86_64 / Windows x86_64（MSVC）  
> 开发语言：Rust stable  
> 参考项目：https://github.com/ulelab/ultraplex  
> 文档日期：2026-08-07

---

## 1. 项目目标

参考 Ultraplex 的核心思想，开发一个**单二进制、跨平台、低依赖、低内存、易部署**的 FASTQ demultiplex 工具，优先解决以下场景：

1. 单端 FASTQ 的 5' barcode demultiplex。
2. paired-end FASTQ 的同步 demultiplex。
3. 5' barcode + 3' barcode 的组合 demultiplex。
4. barcode 中使用 `N` 表示 UMI，并将 UMI 移入 read header。
5. barcode mismatch 容忍。
6. 3' adapter trimming。
7. 5'/3' quality trimming。
8. 最短 read length 过滤。
9. gzip FASTQ 输入与 gzip FASTQ 输出。
10. Linux 和 Windows 使用相同命令和相同行为。

项目重点不是逐行翻译 Ultraplex，而是**重新设计数据流和核心算法**，删除 Python 多进程、临时文件拼接、`cat`/`pigz`/`mv`、SLURM `sbatch` 等平台相关实现。

最终期望：

- `cargo build --release` 后得到一个独立 CLI；
- 普通用户无需 Python、Conda、Cutadapt、pigz 或 SLURM；
- Linux / Windows 行为一致；
- 常规 barcode 数量下内存占用稳定；
- 长 barcode 不再出现 Ultraplex 式 `5^L` reference dictionary 爆炸；
- 结果可通过自动化测试与原版 Ultraplex 进行差异验证。

---

## 2. 对 Ultraplex 的功能理解

Ultraplex 当前主要包含以下处理步骤：

```text
FASTQ input
    ↓
quality trimming
    ↓
3' adapter trimming
    ↓
5' barcode detection
    ↓
5' UMI extraction
    ↓
optional linked 3' barcode detection
    ↓
3' UMI extraction
    ↓
minimum-length filtering
    ↓
sample-specific FASTQ output
```

paired-end 模式还需要：

```text
R1/R2 同步读取
    ↓
利用一端进行 barcode 判断
    ↓
对 mate read 做对应 adapter/barcode trimming
    ↓
R1/R2 必须同时保留或同时丢弃
```

Ultraplex 中值得保留的行为：

- `N` 是 UMI/randomer，不参与 barcode 匹配；
- 5' barcode 的非 `N` 位置必须一致；
- linked 3' barcodes 在同一 5' barcode 分组内位置必须一致；
- barcode mismatch 以最佳匹配为准；
- 如果最佳分数并列，则判定 `no_match`；
- UMI 默认追加到 read header；
- 支持 5' + 3' 组合 barcode；
- 支持只使用 3' barcode 的 paired-end 模式；
- 支持 sample name；
- 支持 `keep barcode`；
- 支持 `ignore/discard no_match`；
- 支持最短 read length 过滤。

不建议照搬的行为：

- 为短 barcode 构建所有可能序列的巨大 reference dictionary；
- Python multiprocessing + Pipe/Queue；
- 每个 worker 输出一批临时文件；
- 最后调用系统 `cat` 合并；
- 外部调用 `pigz`；
- 外部调用 `mv`；
- SLURM `sbatch` compression；
- ultra mode 的超大临时 FASTQ；
- Unix 路径字符串拼接。

---

## 3. 范围冻结

### 3.1 v0.1 必须实现

- FASTQ 4-line record；
- `.fastq` / `.fq`；
- `.fastq.gz` / `.fq.gz`；
- single-end；
- paired-end；
- 5' barcode；
- barcode 中 `N` UMI；
- exact barcode matching；
- barcode mismatch matching；
- linked 3' barcode；
- sample name；
- minimum read length；
- quality trim；
- 3' adapter trim；
- gzip output；
- multithreading；
- Linux；
- Windows；
- summary statistics；
- 与 Ultraplex 的 regression fixtures。

### 3.2 v0.1 明确不实现

以下内容不要在早期 milestone 中加入：

- GUI；
- BAM/SAM；
- FASTA；
- Oxford Nanopore 专用逻辑；
- single-cell barcode correction；
- whitelist learning；
- UMI deduplication；
- interleaved paired FASTQ；
- SLURM；
- `sbatch`；
- 外部 `pigz`；
- Python API；
- network/cloud storage；
- plugin system；
- config language；
- distributed processing。

### 3.3 原 Ultraplex 参数中可以删除的概念

Rust 版直接流式写最终结果，因此以下参数原则上不再需要：

```text
--ultra
--sbatchcompression
--ignore_space_warning
--dont_build_reference
```

其中 `--dont_build_reference` 不再需要，是因为 Rust 版本从设计上就不创建指数级 reference dictionary。

---

## 4. 建议 CLI

开发期二进制暂命名：

```bash
ultraplex-rs
```

基础用法：

```bash
ultraplex-rs demux \
  -i reads.fastq.gz \
  -b barcodes.csv \
  -o results
```

paired-end：

```bash
ultraplex-rs demux \
  -i reads_R1.fastq.gz \
  -I reads_R2.fastq.gz \
  -b barcodes.csv \
  -o results \
  -t 8
```

建议参数：

```text
demux

Required:
  -i, --input <FASTQ>
  -b, --barcodes <CSV>

Paired-end:
  -I, --input2 <FASTQ>

Output:
  -o, --out-dir <DIR>
  -p, --prefix <STR>
      --compression-level <1..9>
      --discard-unassigned

Barcode:
      --mismatches-5 <INT>
      --mismatches-3 <INT>
      --keep-barcodes
      --three-prime-only
      --tso-pattern <PATTERN>

Trimming:
  -a, --adapter-r1 <SEQ>
      --adapter-r2 <SEQ>
  -q, --quality-cutoff-3 <INT>
      --quality-cutoff-5 <INT>
  -l, --min-length <INT>
      --min-adapter-overlap <INT>
      --adapter-error-rate <FLOAT>

Performance:
  -t, --threads <INT>
      --chunk-reads <INT>

Reporting:
      --summary <PATH>
      --quiet
      --log-level <LEVEL>
```

为方便原 Ultraplex 用户迁移，可以在后续加入兼容 alias：

```text
-i2  -> --input2
-m5  -> --mismatches-5
-m3  -> --mismatches-3
-q5  -> --quality-cutoff-5
-kbc -> --keep-barcodes
-inm -> --discard-unassigned
```

不要为了参数完全一致而保留已经没有意义的内部实现参数。

---

## 5. Barcode CSV 兼容策略

优先支持 Ultraplex 当前 CSV 格式。

示例：

```csv
NNNATGNN:sample1,
NNNCCGNN,ATG:sample2,TCA:sample3
NNNCACNN,
```

内部解析后不要继续使用字符串拼接，而是转成结构体。

建议：

```rust
struct BarcodeConfig {
    five_prime: Vec<FivePrimeBarcode>,
    sample_names: HashMap<SampleKey, String>,
}

struct FivePrimeBarcode {
    raw: Vec<u8>,
    informative_positions: Vec<usize>,
    umi_positions: Vec<usize>,
    linked_three_prime: Vec<ThreePrimeBarcode>,
    sample_name: Option<String>,
}

struct ThreePrimeBarcode {
    raw: Vec<u8>,
    informative_positions_from_end: Vec<usize>,
    umi_positions_from_end: Vec<usize>,
    sample_name: Option<String>,
}
```

CSV parser 必须检查：

- 只允许 `A/C/G/T/N`；
- 自动转换为 uppercase；
- sample name 不允许重复；
- barcode 不能为空；
- 5' barcode informative position 一致；
- linked 3' barcode informative position 一致；
- 有 linked 3' barcode 时不能给 5' barcode 分配 sample name；
- mismatch 数不能大于 informative bases 数；
- 输出 filename 中非法字符需要 sanitize；
- Windows 保留文件名如 `CON`, `PRN`, `AUX`, `NUL` 必须处理。

---

## 6. Barcode 匹配算法

### 6.1 不实现 Ultraplex 的 `5^L` reference dictionary

Ultraplex 会枚举：

```text
A C G T N
```

的所有组合。

长度为 `L` 时规模约：

```text
5^L
```

这也是长 barcode 时 reference building 速度和 RAM 急剧恶化的根本原因。

Rust 版禁止使用这一设计。

---

### 6.2 MVP 匹配方式

每个 barcode 编译成：

```rust
CompiledBarcode {
    id,
    raw,
    informative_positions,
    expected_bases,
    umi_positions,
}
```

对 read：

1. 按 informative positions 提取碱基；
2. 与候选 barcode 比较；
3. 计算 Hamming distance；
4. distance <= mismatch threshold 才接受；
5. 找到 distance 最小的 barcode；
6. 最优 distance 出现多个 barcode 时返回 `Unassigned::Ambiguous`。

伪代码：

```text
best_distance = infinity
best_id = none
tie = false

for barcode in candidates:
    d = hamming(observed, barcode.expected)

    if d < best_distance:
        best_distance = d
        best_id = barcode.id
        tie = false
    else if d == best_distance:
        tie = true

if best_distance > allowed_mismatches:
    no_match
else if tie:
    no_match
else:
    best_id
```

barcode 通常很短，因此先保证正确性，再做 SIMD/bit packing。

---

### 6.3 性能优化阶段

MVP 正确后增加：

```text
ExactHashMatcher
PackedHammingMatcher
LinearFallbackMatcher
```

策略：

- mismatch = 0：HashMap exact lookup；
- barcode 长度 <= 32：2-bit packing；
- mismatch > 0：packed Hamming；
- 极长 barcode：普通 byte comparison + early exit。

可将：

```text
A = 00
C = 01
G = 10
T = 11
```

压入 `u64`。

MVP 不要一开始实现复杂 bit trick。

---

## 7. UMI 处理

`N` 不参与 barcode identity matching，而是代表需要提取的 UMI base。

例如：

```text
barcode = NNNATGNN
read    = ACGATGTCxxxxxxxx
```

则：

```text
identity barcode = ATG
UMI              = ACGTC
```

默认 header 兼容 Ultraplex：

```text
@read123rbc:ACGTC
```

建议内部保留：

```rust
struct Umi {
    bases: Vec<u8>,
}
```

若 read header 已存在：

```text
rbc:
```

则沿用原行为，在已有 UMI 后继续追加，而不是重复添加 tag。

同时提供内部测试确保：

- barcode 被 trim 后 sequence/quality 长度一致；
- `--keep-barcodes` 时 sequence 不 trim；
- UMI 仍写入 header；
- 多次处理不会生成多个 `rbc:`。

---

## 8. 5' barcode demultiplex

处理顺序：

```text
read
 ↓
检查 read 长度
 ↓
根据 informative positions 提取 barcode identity
 ↓
barcode matcher
 ↓
提取 N positions 对应 UMI
 ↓
更新 header
 ↓
如果未指定 --keep-barcodes
    trim 整个 barcode pattern 长度
```

注意：

`NNNATGNN` 应删除整个 8 nt，而不是只删除 `ATG`。

---

## 9. 3' barcode demultiplex

3' barcode positions 必须从 read 尾端计算。

例如：

```text
NNATGNNN
```

informative bases 相对于 read 末端的位置必须预先编译，不能每条 read 重算。

建议：

```rust
struct CompiledThreePrimeBarcode {
    pattern_len: usize,
    informative_offsets_from_end: Vec<usize>,
    umi_offsets_from_end: Vec<usize>,
}
```

必须支持：

```text
5' barcode A -> 3' barcode x / y / z
5' barcode B -> 3' barcode p / q
```

因此 3' matcher 应以 5' barcode group 为单位。

---

## 10. Quality trimming

Ultraplex 使用与 BWA `bwa_trim_read` 相同思路的 quality trimming。

Rust 版应直接重写这一小段算法，不需要引入大型依赖。

接口：

```rust
fn quality_trim_bounds(
    qualities: &[u8],
    cutoff_5: u8,
    cutoff_3: u8,
    phred_offset: u8,
) -> (usize, usize)
```

要求：

- 默认 Phred+33；
- 5' cutoff 默认 0；
- 3' cutoff 默认与 CLI 一致；
- sequence 和 quality 使用同一 `(start, end)`；
- 空 read 合法返回；
- 单元测试覆盖原版算法的边界情况。

必须先做 quality trimming，再做 adapter trimming，以保持 Ultraplex 行为接近。

---

## 11. Adapter trimming

不要把 Cutadapt 整个代码库移植到 Rust。

本项目只需要：

```text
3' adapter detection + trimming
```

目标接口：

```rust
struct AdapterMatch {
    read_start: usize,
    read_end: usize,
    adapter_start: usize,
    adapter_end: usize,
    errors: usize,
    overlap: usize,
}

fn find_3p_adapter(
    read: &[u8],
    adapter: &[u8],
    max_error_rate: f32,
    min_overlap: usize,
) -> Option<AdapterMatch>;
```

MVP 算法建议：

1. 优先检测 exact suffix/prefix overlap；
2. 然后对 read 3' 区域运行小规模 semi-global alignment；
3. adapter 通常只有几十 nt，因此 DP 矩阵很小；
4. 支持 substitution；
5. M4 完成前加入 insertion/deletion；
6. 按 overlap 计算允许 error：
   `floor(overlap * max_error_rate)`；
7. 找到最佳合法 match 后 trim read 的 adapter 起点之后内容。

不要在整个 read 上做全局大矩阵。

搜索窗口建议：

```text
read tail length =
adapter_length + max_indel_margin + search_margin
```

后续性能不足时再考虑 Myers bit-parallel matcher。

---

## 12. Single-end 3' barcode 的 `min_trim`

原版逻辑用于避免在 insert 本身的 3' 端偶然匹配 barcode。

Rust 版保留：

```text
只有检测到至少 min_adapter_overlap_for_demux 个 adapter bases，
才允许对 single-end read 使用 3' barcode。
```

建议 CLI：

```text
--min-adapter-trim-for-3p 3
```

不要把它和 adapter matcher 自己的最低 overlap 混成一个参数。

---

## 13. Paired-end 设计

定义：

```rust
struct ReadPair {
    r1: FastqRecord,
    r2: FastqRecord,
}
```

输入层必须验证：

- R1/R2 record 数相同；
- read ID 一致；
- `/1`, `/2` 或 Illumina suffix 允许标准化后匹配；
- 不一致立即报错并给出 record index。

不要静默继续。

输出原则：

```text
如果 pair 不满足 min length：
    R1/R2 都不写

如果 pair 被 sample X 分配：
    R1 -> sample X R1 output
    R2 -> sample X R2 output
```

---

## 14. `--three-prime-only`

保留原 Ultraplex 的核心语义，但不要通过“交换输入文件变量”这种方式实现。

建议显式建模：

```rust
enum DemuxMode {
    FivePrime,
    FiveAndThreePrime,
    ThreePrimeOnly,
}
```

`ThreePrimeOnly` 中：

- 必须 paired-end；
- barcode CSV 以 read1 方向表示；
- matcher 使用 R2；
- 在初始化阶段将 pattern reverse-complement；
- 输出仍然保持真实的 R1/R2 文件身份。

这样可以避免内部变量 R1/R2 互换导致代码难以维护。

---

## 15. TSO pattern

支持与 Ultraplex 相同概念：

```text
NNNNNIII
```

含义：

```text
N -> 加入 UMI
I -> trim 但不加入 UMI
```

解析为：

```rust
struct TsoPattern {
    total_len: usize,
    umi_positions: Vec<usize>,
    ignored_positions: Vec<usize>,
}
```

不要在处理每条 read 时重新分析字符串。

---

## 16. FASTQ I/O

### 16.1 MVP 推荐

优先使用成熟 FASTQ parser，并把其封装在项目内部接口之后。

建议候选：

- `needletail`
- `noodles-fastq`

当前建议：

```text
MVP：needletail
```

原因：

- Rust；
- minimal-copying；
- FASTQ 支持成熟；
- 可处理 gzip；
- Linux/Windows；
- 不需要 Python。

但业务代码绝对不能直接散布 `needletail` 类型。

统一包装：

```rust
trait FastqSource {
    fn next_record(&mut self) -> Result<Option<OwnedFastqRecord>>;
}
```

以后 benchmark 发现 parser 是瓶颈时，可以换后端而不动 demux 逻辑。

### 16.2 record 内部结构

```rust
struct OwnedFastqRecord {
    name: Vec<u8>,
    sequence: Vec<u8>,
    qualities: Vec<u8>,
}
```

不要在热循环中反复转成 UTF-8 `String`。

sequence、quality、barcode 全部优先使用：

```text
&[u8]
Vec<u8>
```

---

## 17. gzip I/O

第一阶段：

```text
flate2
```

默认采用跨平台、易构建方案。

原则：

- 输入自动识别 `.gz`；
- 输出默认 `.fastq.gz`；
- 不调用系统 gzip；
- 不调用 pigz；
- 不调用 shell；
- Windows 与 Linux 共用代码。

后续 benchmark 如果 gzip 是主要瓶颈：

1. benchmark `flate2` 默认 backend；
2. benchmark `zlib-rs` backend；
3. 再决定是否增加 feature：

```toml
[features]
default = ["portable-gzip"]
fast-gzip = []
```

不要在 M0-M5 过早引入 C toolchain 依赖。

---

## 18. 多线程数据流

### 18.1 推荐架构

```text
                 ┌──────────────┐
FASTQ ─────────► │ Reader Thread│
                 └──────┬───────┘
                        │ chunks
                        ▼
              bounded channel
                        │
       ┌────────────────┼────────────────┐
       ▼                ▼                ▼
   Worker 0         Worker 1         Worker N
       │                │                │
       └────────────────┼────────────────┘
                        ▼
               processed chunks
                        │
                        ▼
                 Writer Thread
                        │
                        ▼
              final .fastq.gz files
```

### 18.2 Reader

每次读取固定数量 records：

```text
默认 4096 或 8192 reads/chunk
```

paired-end 则是相同数量 pair。

结构：

```rust
struct InputChunk {
    id: u64,
    records: Vec<ReadOrPair>,
}
```

使用 bounded channel，防止 worker 慢时 reader 无限占内存。

---

### 18.3 Worker

每个 worker：

```text
quality trim
adapter trim
barcode match
UMI extraction
length filter
format FASTQ bytes
local statistics
```

返回：

```rust
struct ProcessedChunk {
    id: u64,
    outputs: HashMap<OutputKey, Vec<u8>>,
    stats: ChunkStats,
}
```

worker 不直接打开输出文件。

这样避免多个线程同时写同一 gzip stream。

---

### 18.4 Writer

Writer 负责：

- 按 sample 打开最终输出；
- 接收 chunk；
- 必要时按 chunk ID reorder；
- 将数据直接写入最终 `.fastq.gz`；
- 关闭并 finish gzip encoder。

不再产生：

```text
_tmp_thread_0.fastq.gz
_tmp_thread_1.fastq.gz
...
```

也不需要最后 concatenate。

---

### 18.5 输出顺序

建议保证：

```text
同一个 sample 内 read 顺序与输入顺序一致
```

实现：

```rust
BTreeMap<chunk_id, ProcessedChunk>
next_expected_chunk
```

writer 只写 `next_expected_chunk`。

这是比“哪个 worker 先完成就先写”更容易测试和复现的设计。

---

## 19. 输出文件设计

建议默认：

single-end：

```text
<prefix>_<sample>.fastq.gz
<prefix>_unassigned.fastq.gz
```

paired-end：

```text
<prefix>_<sample>_R1.fastq.gz
<prefix>_<sample>_R2.fastq.gz
<prefix>_unassigned_R1.fastq.gz
<prefix>_unassigned_R2.fastq.gz
```

没有 sample name：

```text
<prefix>_5bc_<barcode>.fastq.gz
```

组合 barcode：

```text
<prefix>_5bc_<5p>_3bc_<3p>.fastq.gz
```

文件名生成必须集中在：

```text
src/output/naming.rs
```

不要在处理逻辑中拼接路径。

---

## 20. 统计信息

建议：

```rust
struct RunStats {
    total_reads: u64,
    quality_trimmed: u64,
    adapter_trimmed: u64,
    assigned: u64,
    unassigned: u64,
    ambiguous: u64,
    too_short: u64,
    five_prime_matched_three_prime_missing: u64,
    per_sample: HashMap<SampleKey, u64>,
}
```

终端打印：

```text
Total reads:              10,000,000
Assigned:                  9,215,443  (92.15%)
Unassigned:                  612,112   (6.12%)
Ambiguous:                   172,445   (1.72%)
Quality trimmed:           4,312,119  (43.12%)
Adapter trimmed:           8,902,771  (89.03%)
Length filtered:              81,551   (0.82%)
```

同时保存：

```text
<prefix>.summary.tsv
```

后续可选：

```text
<prefix>.summary.json
```

---

## 21. 建议源码结构

保持**单 Cargo package**，不要一开始创建复杂 workspace。

```text
ultraplex-rs/
├── Cargo.toml
├── Cargo.lock
├── LICENSE
├── README.md
├── DEVELOPMENT_PLAN.md
├── CHANGELOG.md
├── src/
│   ├── main.rs
│   ├── lib.rs
│   ├── cli.rs
│   ├── error.rs
│   ├── fastq/
│   │   ├── mod.rs
│   │   ├── reader.rs
│   │   ├── writer.rs
│   │   └── record.rs
│   ├── barcode/
│   │   ├── mod.rs
│   │   ├── config.rs
│   │   ├── matcher.rs
│   │   ├── five_prime.rs
│   │   └── three_prime.rs
│   ├── trim/
│   │   ├── mod.rs
│   │   ├── quality.rs
│   │   └── adapter.rs
│   ├── pipeline/
│   │   ├── mod.rs
│   │   ├── reader.rs
│   │   ├── worker.rs
│   │   └── writer.rs
│   ├── output/
│   │   ├── mod.rs
│   │   └── naming.rs
│   ├── stats.rs
│   └── util.rs
├── tests/
│   ├── fixtures/
│   ├── cli.rs
│   ├── single_end.rs
│   ├── paired_end.rs
│   ├── barcode_matching.rs
│   ├── trimming.rs
│   └── ultraplex_compat.rs
├── benches/
│   ├── barcode.rs
│   ├── adapter.rs
│   └── pipeline.rs
└── .github/
    └── workflows/
        ├── ci.yml
        └── release.yml
```

---

# 22. Milestone 开发计划

---

## M0 — 仓库初始化与跨平台基线

### 目标

只搭建工程，不实现 demultiplex。

### 任务

- 初始化 Cargo；
- Rust edition 使用当前 stable 推荐 edition；
- 添加：
  - `clap`
  - `thiserror`
  - `anyhow`
  - `crossbeam-channel`
  - `csv`
  - `flate2`
  - FASTQ parser；
- 创建上述基础 module；
- CLI 能显示：
  - `--help`
  - `--version`；
- GitHub Actions：
  - Ubuntu；
  - Windows；
- 配置：
  - `cargo fmt --check`
  - `cargo clippy --all-targets --all-features -- -D warnings`
  - `cargo test --all`。

### 禁止

M0 不实现：

- FASTQ processing；
- barcode；
- adapter；
- UMI；
- multithread pipeline。

### 验收

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all
cargo build --release
```

Windows 和 Linux CI 全绿。

---

## M1 — FASTQ 与 gzip I/O

### 目标

可靠读取和写回 FASTQ。

### 任务

实现：

```text
FASTQ reader
FASTQ writer
gzip detection
paired reader
record validation
```

增加命令：

```bash
ultraplex-rs validate -i reads.fastq.gz
```

输出：

```text
records
bases
min length
max length
paired consistency
```

### 测试

- uncompressed；
- gzip；
- empty file；
- truncated record；
- sequence/quality 长度不同；
- paired count 不一致；
- paired ID 不一致；
- Windows CRLF；
- spaces in read name。

### 验收

对 fixture：

```text
read -> write -> read
```

record 内容完全一致。

---

## M2 — Barcode CSV + 5' exact demultiplex

### 目标

完成最小可用 demultiplexer。

仅实现：

```text
single-end
5' barcode
mismatch = 0
UMI
sample name
keep barcode
unassigned
```

### 任务

- CSV parser；
- barcode validation；
- compiled barcode；
- exact matcher；
- 5' UMI extraction；
- output naming；
- sample writer；
- basic stats。

### 验收 fixture

barcode：

```text
NNNATGNN:sampleA
NNNCCGNN:sampleB
```

确保：

- sampleA 正确；
- sampleB 正确；
- no_match 正确；
- UMI 正确；
- trim 长度正确；
- quality 同步 trim；
- output record 顺序不变。

---

## M3 — Barcode mismatch + linked 3' demultiplex

### 目标

完成 Ultraplex 最核心特色。

### 任务

- Hamming matcher；
- mismatch threshold；
- ambiguous tie；
- linked barcode groups；
- 3' position extraction；
- 3' UMI；
- combined `SampleKey`；
- `--mismatches-5`；
- `--mismatches-3`。

### 必测

```text
perfect match
1 mismatch
beyond threshold
two equally good barcodes
read contains N
short read
3p linked match
5p match but 3p no match
```

### 验收

所有 barcode 决策必须是 deterministic。

---

## M4 — Quality trim + adapter trim

### M4.1 Quality trim

按 Ultraplex/BWA 算法实现。

必须先写 unit tests，然后接 pipeline。

### M4.2 Exact adapter trim

先实现：

```text
exact 3' adapter
partial adapter
minimum overlap
```

### M4.3 Approximate adapter trim

加入：

```text
substitution
insertion
deletion
max error rate
```

实现 small semi-global DP。

### M4.4 Single-end 3' gating

实现：

```text
min adapter trim -> 是否允许 3' barcode demux
```

### 验收

建立 20–50 个固定 adapter trimming fixture，对照 Ultraplex 输出。

---

## M5 — Paired-end + three-prime-only + TSO

### 目标

完成主要功能 parity。

### 任务

- paired synchronized reader；
- paired filtering；
- R1/R2 output；
- mate trimming；
- `--three-prime-only`；
- reverse-complement barcode；
- TSO `N/I` pattern；
- paired stats。

### 验收

至少建立：

```text
normal paired
R1/R2 ID mismatch
R1/R2 count mismatch
5p only
5p + 3p
3p only
TSO
short R1
short R2
```

---

## M6 — Chunked multithreading

### 目标

在不改变结果的前提下并行化。

### 任务

实现：

```text
Reader -> workers -> ordered writer
```

使用：

```text
bounded crossbeam channels
```

要求：

- `--threads 1` 与 `--threads 8` 输出解压后完全一致；
- sample 内顺序一致；
- 无临时文件；
- Ctrl+C/worker error 不留下未 finish 的假正常 gzip；
- pipeline 任何线程 error 都传回主线程。

### 内存目标

默认：

```text
< 512 MB
```

对于普通 Illumina read + 8 threads + 默认 chunk。

内存不应随输入 FASTQ 总大小增长。

---

## M7 — Differential testing against Ultraplex

这是整个项目最重要的质量阶段之一。

### 方法

编写 synthetic data generator：

输入：

```text
barcode config
read length
error rate
adapter presence
UMI
mismatch positions
pair mode
```

生成：

```text
R1.fastq.gz
R2.fastq.gz
barcodes.csv
expected metadata.tsv
```

同时运行：

```bash
ultraplex ...
ultraplex-rs ...
```

比较：

```text
sample assignment
UMI
sequence
quality
trim positions
read counts
unassigned
```

gzip 二进制本身无需一致。

比较方式：

```text
gunzip -> normalized FASTQ -> diff
```

允许明确记录的兼容差异，但必须写入：

```text
docs/COMPATIBILITY.md
```

禁止出现“测试不一致但先忽略”的情况。

---

## M8 — 性能优化

先 profile，后优化。

### Benchmark dataset

至少准备：

```text
1 million reads
10 million reads
paired 10 million reads
4 barcodes
32 barcodes
96 barcodes
384 barcodes
barcode length 6
barcode length 10
barcode length 16
```

线程：

```text
1
2
4
8
16
```

记录：

```text
reads/sec
MB/sec compressed input
wall time
peak RSS
output size
CPU utilization
```

Linux 与 Windows 都运行。

### 优化顺序

只按 profile 结果行动：

1. allocations；
2. FASTQ parsing；
3. barcode matching；
4. adapter alignment；
5. gzip compression；
6. output HashMap；
7. chunk size。

### 可能优化

- reuse buffers；
- byte slice 替代 String；
- exact barcode HashMap；
- 2-bit packed barcode；
- early-exit Hamming；
- preallocated sample buffers；
- `SmallVec`；
- faster gzip backend feature；
- custom strict FASTQ parser；
- Myers adapter matching。

不要未经 benchmark 就使用 SIMD/unsafe。

---

## M9 — Release

### GitHub Actions release matrix

至少：

```text
x86_64-pc-windows-msvc
x86_64-unknown-linux-gnu
```

可选：

```text
x86_64-unknown-linux-musl
aarch64-unknown-linux-gnu
```

release 包：

Windows：

```text
ultraplex-rs-vX.Y.Z-x86_64-pc-windows-msvc.zip
```

Linux：

```text
ultraplex-rs-vX.Y.Z-x86_64-unknown-linux-gnu.tar.gz
```

每个 release 包包含：

```text
binary
README.md
LICENSE
CHANGELOG.md
```

---

# 23. 测试策略

测试比例建议：

```text
Unit tests          50%
Integration tests   30%
Compatibility tests 15%
Benchmarks           5%
```

核心 module 都必须有 unit tests。

---

## 23.1 Property tests

建议加入 `proptest`，尤其用于：

### barcode

性质：

```text
distance(a, a) = 0
distance(a, b) = distance(b, a)
exact barcode 必须命中自己
超过 threshold 必须 no_match
tie 必须 ambiguous
```

### quality trim

性质：

```text
0 <= start <= end <= len
seq.len == qual.len
```

### reverse complement

性质：

```text
revcomp(revcomp(x)) == x
```

---

## 23.2 Fuzzing

v0.1 之后加入：

```text
cargo-fuzz
```

targets：

```text
FASTQ parser
barcode CSV parser
adapter matcher
```

重点防：

- panic；
- integer underflow；
- malformed UTF-8；
- truncated FASTQ；
- extreme barcode length。

---

# 24. Windows 特别要求

所有开发 milestone 都必须同时考虑 Windows。

禁止：

```rust
Command::new("cat")
Command::new("mv")
Command::new("pigz")
```

禁止：

```text
"/" 手动拼路径
Unix file descriptor 假设
fork()
```

统一使用：

```rust
Path
PathBuf
std::fs
std::io
```

CI 必须使用：

```text
windows-latest
```

测试 output path：

```text
C:\path with spaces\test data\
```

还要测试：

- Unicode directory；
- long-ish path；
- CRLF CSV；
- Windows filename reserved names；
- existing output directory；
- overwrite policy。

---

# 25. Error handling

不要使用大量：

```rust
unwrap()
expect()
panic!()
```

处理用户输入。

定义：

```rust
enum AppError {
    Io,
    FastqFormat,
    BarcodeConfig,
    InvalidPair,
    InvalidBase,
    OutputConflict,
    WorkerFailure,
}
```

CLI error 示例：

```text
error: paired FASTQ records are out of sync
  R1 record 18291: A01234:87:H...
  R2 record 18291: A01234:91:H...
```

而不是：

```text
thread panicked at ...
```

---

# 26. 输出覆盖策略

默认：

```text
如果目标文件存在 -> 报错退出
```

增加：

```text
--force
```

才允许覆盖。

禁止像原程序一样启动后自动清理与 prefix 匹配的旧文件。

这是数据安全要求。

---

# 27. 日志

普通状态写 stderr。

默认：

```text
INFO
```

参数：

```text
--quiet
--log-level warn|info|debug|trace
```

不要让 worker 直接大量打印。

每隔固定 read 数或固定时间由主线程输出：

```text
processed 10.0 M reads | 1.82 M reads/s
```

---

# 28. 性能目标

这些是工程目标，不是预先承诺。

### Correctness 优先

首先：

```text
Rust output == expected output
```

然后再优化速度。

### v0.1 目标

对常见 Illumina demultiplex：

```text
8 threads:
>= 原 Ultraplex 速度的 70%
```

即可接受首版。

### v0.2 目标

经过 profiling：

```text
达到或超过 Ultraplex
```

### 内存目标

barcode 数量常规时：

```text
< 512 MB
```

并且：

```text
RAM ≈ chunk_size × in-flight chunks
```

而不是与 FASTQ 总大小或 `5^barcode_length` 成正比。

---

# 29. 依赖策略

优先：

```text
small
mature
cross-platform
actively maintained
```

建议初始依赖：

```text
clap
anyhow
thiserror
crossbeam-channel
csv
flate2
needletail
serde (仅需要结构化 report 时)
```

开发依赖：

```text
tempfile
assert_cmd
predicates
proptest
criterion
```

原则：

- 不引入 Tokio；
- 不需要 async；
- 不引入完整 `rust-bio` 作为 MVP 必需依赖；
- 不引入 C/C++ dependency 作为默认路径；
- Cargo.lock 提交到仓库；
- release binary 尽可能自包含。

如果后续 benchmark 证明 adapter matcher 是瓶颈，可评估 Rust-Bio Myers 或内部 bit-parallel 实现。

---

# 30. License

Ultraplex 使用 MIT License，并且仓库明确说明其中大量代码来源于 Cutadapt，后者相应代码同样附带 MIT 许可。

建议新项目：

```text
MIT
```

如果只是根据行为和算法重新实现，保留清晰 reference 即可。

如果直接复制、翻译或修改 Ultraplex/Cutadapt 的实质性源码片段，则必须按原许可证要求保留对应 copyright 和 permission notice。

建议增加：

```text
NOTICE.md
```

写明：

```text
This project was inspired by Ultraplex:
https://github.com/ulelab/ultraplex

Ultraplex contains code derived from Cutadapt:
https://github.com/marcelm/cutadapt
```

---

# 31. Codex / Grok Build 执行规则

将本文件交给 coding agent 时，使用以下规则。

## 总原则

```text
一次只执行一个 Milestone。
```

Agent 不得自行跳到下一个 milestone。

每个 milestone 完成必须：

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all
cargo build --release
```

全部通过后才能结束该 milestone。

---

## 每个 milestone 的提交要求

建议 commit：

```text
M0: initialize cross-platform Rust CLI
M1: add FASTQ and gzip IO
M2: implement 5-prime barcode demultiplexing
M3: add mismatch and linked 3-prime barcodes
M4: add quality and adapter trimming
M5: add paired-end and three-prime-only modes
M6: add chunked multithread pipeline
M7: add Ultraplex differential tests
M8: optimize hot paths
M9: add release packaging
```

禁止一个 commit 横跨多个 milestone。

---

# 32. 推荐给 Coding Agent 的第一条任务

以下内容可直接交给 Codex / Grok Build：

```text
执行 DEVELOPMENT_PLAN.md 中的 M0。

严格范围冻结：
只完成 Rust 仓库初始化、模块骨架、CLI --help/--version、
Linux + Windows CI、fmt/clippy/test/build 基线。

不要实现 FASTQ parser。
不要实现 barcode。
不要实现 UMI。
不要实现 adapter trimming。
不要实现 quality trimming。
不要实现 multithreading pipeline。
不要提前执行 M1 或之后的任务。

完成后运行：

cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all
cargo build --release

如果任何命令失败，修复到全部通过。

最终报告：
1. 创建/修改了哪些文件；
2. dependency 列表；
3. Linux/Windows CI 配置；
4. 四条验证命令的结果；
5. 明确声明没有实现 M1+ 功能。
```

---

# 33. M0 完成后的第二条任务

```text
执行 DEVELOPMENT_PLAN.md 中的 M1。

只实现 FASTQ/gzip I/O 与 validate 命令。

范围：
- FASTQ read/write
- gzip read/write
- paired record synchronization
- malformed FASTQ errors
- validate command
- fixtures and tests

不要实现：
- barcode
- UMI
- adapter
- quality trimming
- demultiplex
- multithread processing

要求 Linux/Windows 均通过：

cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all
cargo build --release
```

---

# 34. 最终完成标准

项目只有同时满足以下条件才算达到 v0.1：

- [ ] Windows release binary 可直接运行
- [ ] Linux release binary 可直接运行
- [ ] single-end demux
- [ ] paired-end demux
- [ ] gzip input
- [ ] gzip output
- [ ] 5' barcode
- [ ] 3' barcode
- [ ] linked combinatorial barcode
- [ ] UMI extraction
- [ ] mismatch handling
- [ ] ambiguous barcode handling
- [ ] quality trimming
- [ ] adapter trimming
- [ ] min length
- [ ] sample names
- [ ] three-prime-only
- [ ] TSO pattern
- [ ] no external shell commands
- [ ] no temporary per-thread FASTQ architecture
- [ ] bounded-memory multithreading
- [ ] deterministic output
- [ ] differential fixtures against Ultraplex
- [ ] Windows CI
- [ ] Linux CI
- [ ] release artifacts
- [ ] README usage documentation
- [ ] LICENSE / NOTICE complete

---

# 35. 推荐开发优先级

实际 coding 顺序必须保持：

```text
Correct FASTQ I/O
        ↓
Correct barcode model
        ↓
Correct 5' demux
        ↓
Correct 3' demux
        ↓
Correct trimming
        ↓
Correct paired-end
        ↓
Compatibility tests
        ↓
Multithreading
        ↓
Profiling
        ↓
Optimization
        ↓
Release
```

不要采用：

```text
先追求“400M reads 很快”
然后再补 correctness
```

因为 demultiplex 工具最危险的问题不是崩溃，而是**快速地产生看起来正常、实际上分错 sample 的 FASTQ**。

---

# 36. 关键架构决定摘要

最终推荐架构：

```text
Rust stable
single Cargo package
byte-oriented FASTQ records
streaming gzip
compiled barcode patterns
Hamming barcode matching
small semi-global adapter matcher
chunked bounded pipeline
N worker threads
1 ordered writer
direct final output
no temporary FASTQ
no shell commands
no pigz dependency
no SLURM dependency
Linux + Windows CI from M0
```

这是相对于直接“Rust 翻译 Ultraplex”更适合长期维护的路线。

---

## 参考资料

1. Ultraplex  
   https://github.com/ulelab/ultraplex

2. Ultraplex README  
   https://github.com/ulelab/ultraplex/blob/master/README.md

3. Cutadapt  
   https://github.com/marcelm/cutadapt

4. needletail  
   https://docs.rs/needletail/

5. flate2  
   https://docs.rs/flate2/

6. noodles-fastq（备用 FASTQ backend）  
   https://docs.rs/noodles-fastq/

7. Rust-Bio Myers approximate matching（后续优化候选）  
   https://docs.rs/bio/latest/bio/pattern_matching/myers/

---

## 结论

不建议把 Ultraplex 的 Python 多进程和临时文件架构机械翻译成 Rust。

最值得继承的是：

```text
barcode/UMI 语义
5' + 3' combinatorial demultiplex
paired-end workflow
quality + adapter trimming
```

最值得重新设计的是：

```text
barcode matching
FASTQ I/O
threading
compression
output
cross-platform path handling
error handling
testing
```

如果按 M0 → M9 顺序严格推进，这个项目可以先得到一个非常小但正确的 MVP，再逐步达到 Ultraplex 的主要功能，并利用 Rust 在内存安全、线程模型、部署和 Windows 支持上的优势形成一个更轻量的长期工具。
