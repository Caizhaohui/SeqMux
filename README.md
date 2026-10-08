# SeqMux

**SeqMux** is a high-performance, single-binary FASTQ demultiplexer written in Rust.

It demultiplexes single-end and paired-end FASTQ using a **sample barcode table** (dual or single barcode), featuring integrated 3′ adapter and quality trimming, automatic handling of mixed insert orientations in amplicon sequencing, and ordered multi-threaded streaming I/O. No Python, Conda, pigz, or external tools required.

---

## What's New in v0.4

- 🚀 **High-Throughput `mimalloc` Build (~1.65× to 2.1× Speedup)**: Optional `--features mimalloc-allocator` drastically reduces memory allocator lock contention in multi-threaded gzip pipelines, reaching **>720,000 pairs/s** on full-scale sequencing runs with bounded memory footprint (~466 MB peak RSS on 140M+ pairs).
- 🛡️ **Global Barcode Confidence Margin (`--min-mismatch-delta`)**: Enforces a strict separation margin against the closest distinct candidate sample across the entire barcode pool (independent of mismatch thresholds), guarding against ambiguous assignments while collapsing canonical/swapped orientations.
- ⚡ **Zero-Overhead Hot Path**: Preserves the $O(1)$ `ExactDual8Matcher` fast path for standard exact-match production workloads (`mismatches=0, min-mismatch-delta=0`).
- 📊 **Refined Ambiguity Accounting**: Clarifies that ambiguous reads are tracked as an explicit subset of unassigned reads (`total = assigned + unassigned`, with `ambiguous <= unassigned`).
- ✅ **100% Decompressed Output Parity**: Fully audited on real multi-million-pair sequencing runs—standard and mimalloc builds generate byte-for-byte identical decompressed FASTQ data across all samples.

---

## Performance & Comparison with Peer Tools

### 1. Feature & Throughput Comparison

SeqMux was benchmarked against peer demultiplexers (**Ultraplex**, **fqtk**, and **Cutadapt**) on real paired-end amplicon sequencing datasets using a 16-CPU allocation (`x86_64` Linux).

| Tool | Core Language | Throughput (pairs/s) | Peak RSS | Mixed PE Orientations | Integrated Trimming | Barcode Confidence Margin | Single Binary |
|:---|:---:|---:|---:|:---:|:---:|:---:|:---:|
| **SeqMux v0.4 (mimalloc)** | **Rust** | **570,000 – 721,000** | **~466 MB** | **Native auto-canonicalize** | **Yes (adapter + quality)** | **Yes (`--min-mismatch-delta`)** | **Yes** |
| **SeqMux v0.4 (standard)** | **Rust** | **346,000** | **~115 – 155 MB** | **Native auto-canonicalize** | **Yes (adapter + quality)** | **Yes (`--min-mismatch-delta`)** | **Yes** |
| **fqtk (v0.4.1)** | Rust | ~460,000 – 474,000 | ~1,280 MB | Fixed structure only | No (requires external step) | Yes (read structure) | Yes |
| **Ultraplex (v1.2)** | Python | ~24,000 – 29,000 | ~58 MB | Workaround (doubled CSV + merge) | Yes (adapter + quality) | No | No (Python env) |
| **Cutadapt** | Python/C | ~15,000 – 25,000 | ~60 MB | Manual multi-step scripts | Yes (adapter + quality) | No | No (Python env) |

#### Key Advantages:
- **vs. Ultraplex**: SeqMux is **10× to 36× faster**. Ultraplex requires minutes pre-generating $5^L$ reference dictionaries for combinatorial 8-mer barcodes and cannot natively demultiplex mixed PE insert orientations without doubling the barcode table and manually concatenating output FASTQs.
- **vs. fqtk**: In high-throughput multi-threaded pipelines, SeqMux with `mimalloc-allocator` achieves **~1.5× higher throughput** (721k vs 474k pairs/s) while using **~63% less memory** (466 MB vs 1,280 MB peak RSS). Furthermore, SeqMux natively handles mixed insert orientations and 3′ adapter trimming in a single streaming pass, whereas fqtk only supports fixed-orientation read structures without trimming.
- **vs. Cutadapt / In-House Python Scripts**: Replaces complex, error-prone shell scripts with a single portable binary, processing 140+ million pairs in minutes instead of hours.

### 2. Full-Dataset Production Benchmark

Measured on 16-CPU server allocations (`-t 12`, `--compression-level 1`, exact matching, Illumina adapter trimming enabled, both orientations enabled):

| Dataset | Read Pairs | Wall Clock (Standard) | Wall Clock (v0.4 mimalloc) | Speedup | Peak RSS | Assigned Pairs | Result Parity |
|:---|---:|---:|---:|---:|---:|---:|:---:|
| **Library A (35-plex dual-index)** | 140,771,720 | 6 min 47 s (406.8 s) | **3 min 15 s (195.2 s)** | **2.08×** | 466 MB | 87,477,095 (62.1%) | 100% exact match |
| **Library B (18-plex dual-index)** | 60,523,000 | 2 min 57 s (177.0 s) | **1 min 26 s (86.1 s)** | **2.06×** | 430 MB | 39,455,648 (65.2%) | 100% exact match |

*Note: Assignment counts match independent exact demultiplexing references across all samples with 100% parity. Decompressed FASTQ outputs between standard and mimalloc builds are byte-for-byte identical.*

---

## Installation

### Pre-built Binaries
Download pre-compiled binaries for Linux and Windows from the [GitHub Releases](https://github.com/Caizhaohui/SeqMux/releases) page.

### Build from Source

Requirements: Rust stable (1.74+).

```bash
git clone https://github.com/Caizhaohui/SeqMux.git
cd SeqMux

# Standard release build:
cargo build --release

# High-throughput release build with mimalloc (~1.65x - 2.1x speedup):
cargo build --release --features mimalloc-allocator

# The compiled binary is at target/release/seqmux
```

---

## Quick Start

**Paired-end (recommended for dual-barcode amplicon sequencing):**

```bash
seqmux demux \
  -i sample_R1.fq.gz \
  -I sample_R2.fq.gz \
  -b sample_barcodes.csv \
  -o results \
  -p seqmux \
  -t 12 \
  --compression-level 1
```

**Single-end:**

```bash
seqmux demux \
  -i reads.fastq.gz \
  -b sample_barcodes.csv \
  -o results \
  -t 8
```

**Validate FASTQ integrity and pairing:**

```bash
seqmux validate -i sample_R1.fq.gz -I sample_R2.fq.gz
```

---

## Sample Barcode Table (Required Format)

SeqMux reads a **CSV file with a header row**. This is the only supported barcode input format.

### Minimal Example

```csv
SampleNumber,Barcode1,Barcode2
Sample_01,AAGTCCAA,GGAGTACT
Sample_02,GACCTGAA,GGACTTGG
Sample_03,GGCTTAAG,TGGATCGA
```

### Full Example (Extra Columns are Ignored)

Standard laboratory sample sheets with auxiliary columns are supported directly without pre-filtering:

```csv
SampleNumber,Barcode1,Barcode2,PCR_product,rawdata1,rawdata2,library_round,F_primer,R_primer
Sample_01,AAGTCCAA,GGAGTACT,amplicon_1,sample_R1.fq.gz,sample_R2.fq.gz,A,fwd...,rev...
Sample_02,GACCTGAA,GGACTTGG,amplicon_2,sample_R1.fq.gz,sample_R2.fq.gz,A,fwd...,rev...
```

An example template is available at [`examples/sample_barcodes.csv`](examples/sample_barcodes.csv).

### Required Columns & Aliases

| Logical Field | Required | Accepted Column Headers (Case-Insensitive) |
|:---|:---:|:---|
| **Sample ID** | **Yes** | `SampleNumber`, `Sample`, `Sample_Name`, `Sample_ID`, `Name` |
| **Barcode 1** | **Yes** | `Barcode1`, `Barcode_1`, `BC1`, `Index1`, `i7` |
| **Barcode 2** | Optional* | `Barcode2`, `Barcode_2`, `BC2`, `Index2`, `i5` |

\* *If `Barcode2` is present in any row, all rows must specify `Barcode2`. If omitted entirely, SeqMux automatically operates in **single-barcode** mode.*

### Table Rules
- Header row is **required**.
- Sample names must be unique.
- Barcode combinations `(Barcode1, Barcode2)` must be unique.
- Valid bases: `A/C/G/T/N` (case-insensitive). `N` represents a UMI position (extracted to the read header tag `rbc:`).

---

## Barcode Matching & Trimming Pipeline

### Matching Modes

| Mode | Trigger | Matching Logic |
|:---|:---|:---|
| **Dual Barcode + PE** | `Barcode2` present and `-I` provided | **Both orientations (default)**: `Barcode1`@R1 5′ + `Barcode2`@R2 5′, **and** `Barcode2`@R1 5′ + `Barcode1`@R2 5′. Swapped reads are automatically rotated so output R1 carries Barcode1. |
| **Dual Barcode + SE** | `Barcode2` present, no `-I` | `Barcode1` @ R1 5′, `Barcode2` @ R1 3′ |
| **Single Barcode** | No `Barcode2` column | `Barcode1` @ R1 5′ only |

### Processing Order
1. **Barcode Matching**: Evaluates the original 5′ sequences before any trimming.
2. **Barcode Trimming**: Strips matched barcode bases from the reads (unless `--keep-barcodes`).
3. **Quality Trimming**: BWA-style Phred quality trimming from the 3′ end (optional, default `-q 0`).
4. **Adapter Trimming**: Integrated 3′ Illumina adapter trimming (on by default: R1 `AGATCGGAAGAGCACACGTCTGAA`, R2 `AGATCGGAAGAGCGTCGTG`; disable with `--no-adapter`).
5. **Length Filtering**: Discards pairs where either mate falls below `--min-length` (default 0).

### Confidence Margin (`--min-mismatch-delta`)

When `--min-mismatch-delta <INT>` is set ($> 0$):
- The best eligible candidate must satisfy the per-end thresholds (`d1 <= max1`, `d2 <= max2`).
- SeqMux checks the distance to the **closest distinct competitor sample** across the complete sample pool (which does not need to satisfy eligibility).
- For paired-end data, canonical and swapped orientations within the same sample are collapsed before inter-sample comparison, preventing a sample from competing with itself.
- Reads failing the confidence margin are classified into the `ambiguous` subset and routed to unassigned output.
- Setting `--min-mismatch-delta 0` (default) preserves 100% backward compatibility with exact tie-breaking.

---

## Common CLI Options

| Option | Description | Default |
|:---|:---|:---:|
| `-i / --input` | Forward / R1 / single-end FASTQ (`.fq`, `.fastq`, `.gz`) | *Required* |
| `-I / --input2` | Reverse / R2 FASTQ for paired-end sequencing | *Optional* |
| `-b / --barcodes` | Sample barcode table CSV file | *Required* |
| `-o / --out-dir` | Output directory | *Required* |
| `-p / --prefix` | Output file prefix | `seqmux` |
| `-t / --threads` | Number of worker threads | `4` |
| `--compression-level` | Gzip compression level (`1`–`9`, `1` recommended for speed) | `6` |
| `--mismatches-1` | Maximum allowed mismatches for Barcode 1 | `0` |
| `--mismatches-2` | Maximum allowed mismatches for Barcode 2 | `0` |
| `--min-mismatch-delta` | Minimum confidence distance margin to closest runner-up sample | `0` |
| `--orientation` | Dual-barcode PE orientation: `both`, `canonical`, `swapped` | `both` |
| `--no-canonicalize` | Do not rotate swapped mates to Barcode1-on-R1 | `false` |
| `--keep-barcodes` | Keep barcode bases on reads (do not trim) | `false` |
| `--discard-unassigned` | Do not write unassigned read files | `false` |
| `-a / --adapter-r1` | 3′ adapter sequence for R1 | Illumina standard |
| `--adapter-r2` | 3′ adapter sequence for R2 | Illumina standard |
| `--no-adapter` | Disable 3′ adapter trimming | `false` |
| `-q / --quality-cutoff-3` | 3′ Phred quality cutoff | `0` |
| `-l / --min-length` | Minimum read length after trimming | `0` |
| `--counts-only` | Count assignments only; do not write FASTQ files | `false` |
| `--force` | Overwrite existing files in output directory | `false` |

Run `seqmux demux --help` for the full parameter reference.

---

## Output Layout

### Single-end:
```text
<prefix>_<SampleNumber>.fastq.gz
<prefix>_unassigned.fastq.gz
<prefix>.summary.tsv
```

### Paired-end:
```text
<prefix>_<SampleNumber>_R1.fastq.gz
<prefix>_<SampleNumber>_R2.fastq.gz
<prefix>_unassigned_R1.fastq.gz
<prefix>_unassigned_R2.fastq.gz
<prefix>.summary.tsv
```

*Example*: Sample `Sample_01` $\to$ `seqmux_Sample_01_R1.fastq.gz` and `seqmux_Sample_01_R2.fastq.gz`.

---

## Design Principles

- **Streaming Ordered Writer**: Reader $\to$ bounded chunks $\to$ worker threads $\to$ single ordered writer. Worker threads never open output files; the writer handles gzip compression sequentially, ensuring deterministic output record order without temporary intermediate files.
- **Bounded Memory**: Streaming architecture guarantees that memory consumption remains strictly constant regardless of input file size (hundreds of megabytes, not gigabytes).
- **Correctness First**: Demultiplexing prioritizes precision over speed; ambiguous reads are never silently misassigned.

---

## Development & Verification

```bash
# Code formatting check
cargo fmt --check

# Clippy linter check
cargo clippy --all-targets --all-features -- -D warnings

# Run all unit, integration, and property tests
cargo test --all

# Build release binaries
cargo build --release
cargo build --release --features mimalloc-allocator
```

---

## License

MIT License — see [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md).

## Acknowledgements

Pipeline and trimming concepts were informed by demultiplexing practices in the NGS bioinformatics community, including Ultraplex.
