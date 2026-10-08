# SeqMux

**Fast, single-binary FASTQ demultiplexer in Rust for amplicon sequencing.**

SeqMux assigns single-end or paired-end FASTQ reads to samples based on a sample sheet CSV with headers. It supports inline barcodes, mixed paired-end insert orientations, barcode trimming, UMI extraction, Phred quality trimming, 3′ adapter trimming, and length filtering. The tool operates using streaming processing and an ordered writer; running demultiplexing does not require Python, Conda, pigz, or external demultiplexing utilities.

> [!IMPORTANT]
> **Current Scope & Application Notice**:
> SeqMux is specifically engineered and rigorously production-validated for **amplicon sequencing data** (e.g., multiplexed PCR amplicon libraries with inline 5′/3′ barcodes and ~50/50 mixed insert orientations). It has **not yet been benchmarked or tested on other sequencing library types** (such as WGS, RNA-seq, ChIP-seq, or single-cell/spatial platforms). Validation and presets for broader library architectures are planned for future development.

## Table of Contents

- [Applicable Scenarios](#applicable-scenarios)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Sample Table](#sample-table)
- [Matching and Orientations](#matching-and-orientations)
- [Processing Pipeline Order](#processing-pipeline-order)
- [UMI and TSO](#umi-and-tso)
- [Common Usage Examples](#common-usage-examples)
- [Full Parameter Reference](#full-parameter-reference)
- [Outputs and Statistics](#outputs-and-statistics)
- [Performance and Resource Usage](#performance-and-resource-usage)
- [Verification and Development](#verification-and-development)
- [Known Limitations](#known-limitations)
- [Frequently Asked Questions (FAQ)](#frequently-asked-questions-faq)
- [License and Acknowledgements](#license-and-acknowledgements)

---

## Applicable Scenarios

The primary validated use case is multiplexed PCR amplicon FASTQ demultiplexing with inline barcodes, especially paired-end libraries where Barcode1 and Barcode2 may appear on R1/R2 in either orientation or swapped orientations:

| Data Structure | Supported Mode |
|:---|:---|
| PE: Barcode1 at R1 5′, Barcode2 at R2 5′ | Dual barcode, `canonical` orientation |
| PE: Barcode2 at R1 5′, Barcode1 at R2 5′ | Dual barcode, `swapped` orientation |
| Mixed orientations within the same PE library | Default `--orientation both` |
| SE: Barcode1 at 5′, Barcode2 at 3′ | Dual barcode SE mode |
| SE or PE: only Barcode1 at R1 5′ | Single barcode mode |
| Barcode with embedded UMI | Mark UMI positions with `N` in CSV |

> [!NOTE]
> `Barcode1` and `Barcode2` refer to sequences located directly within the sequenced reads. Do not pass sequencer sample sheet i7/i5 index tables as inline barcode tables. While column aliases accept `i7`/`i5`, this does not mean the tool automatically reads separate index FASTQ files or performs i5 reverse-complement conversion.
> 
> WGS, RNA-seq, ChIP-seq, single-cell, spatial transcriptomics, and other read structures have not received equivalent production validation. The tool does not search for barcodes at arbitrary positions, does not perform indel-tolerant edit distance error correction on barcodes, and is not a general-purpose read-structure parser.

---

## Installation

### Build from Source

Requires Rust 1.85 or higher and platform build tools. Linux and Windows are covered by CI.

```bash
git clone https://github.com/Caizhaohui/SeqMux.git
cd SeqMux
cargo build --release --locked
./target/release/seqmux --version
```

On Windows, the executable is generated at `target/release/seqmux.exe`. Multi-line commands below use Bash continuation syntax (`\`); in PowerShell, enter commands as single lines or use the backtick (`` ` ``) continuation character.

### Optional `mimalloc` Build

```bash
cargo build --release --locked --features mimalloc-allocator
```

`mimalloc-allocator` is a cargo compile-time feature, not a runtime flag of `seqmux demux`. It provides throughput gains on Linux HPC workloads but may increase peak RSS. The default build uses the system allocator.

To retain both binaries for benchmarking, use separate target directories:

```bash
cargo build --release --locked --target-dir target/standard
cargo build --release --locked --features mimalloc-allocator --target-dir target/mimalloc

./target/standard/release/seqmux --version
./target/mimalloc/release/seqmux --version
```

### Precompiled Packages

Download packages for your operating system from [GitHub Releases](https://github.com/Caizhaohui/SeqMux/releases). Current release workflows target `x86_64` Linux GNU and `x86_64` Windows MSVC using the default allocator. If mimalloc is desired, build from source as shown above. Verify the actual version with `--version`.

---

## Quick Start

Prepare `barcodes.csv`:

```csv
SampleNumber,Barcode1,Barcode2
Sample_01,AAGTCCAA,GGAGTACT
Sample_02,GACCTGAA,GGACTTGG
Sample_03,GGCTTAAG,TGGATCGA
```

Validate paired-end inputs:

```bash
seqmux validate -i reads_R1.fastq.gz -I reads_R2.fastq.gz
```

Run demultiplexing:

```bash
seqmux demux \
  -i reads_R1.fastq.gz \
  -I reads_R2.fastq.gz \
  -b barcodes.csv \
  -o results \
  -p run01 \
  -t 12 \
  --compression-level 1
```

This command runs exact barcode matching, evaluates both canonical and swapped orientations, canonicalizes swapped matches so output R1 carries Barcode1, strips matched barcodes, and applies default 3′ adapter trimming. Quality trimming defaults to off; `--min-length` defaults to 0.

`-t 12` is an empirical setting from HPC runs; `--compression-level 1` explicitly selects fast compression (default is 6 when omitted).

---

## Sample Table

### Format and Column Names

The only native sample sheet format is CSV with a header row. Auxiliary columns in laboratory sample sheets are ignored.

| Logical Field | Required | Common Accepted Column Names (Case-Insensitive) |
|:---|:---:|:---|
| Sample ID | Required | `SampleNumber`, `Sample`, `Sample_Name`, `Sample_ID`, `Name` |
| Barcode1 | Required | `Barcode1`, `Barcode_1`, `BC1`, `Index1`, `i7` |
| Barcode2 | Optional | `Barcode2`, `Barcode_2`, `BC2`, `Index2`, `i5` |

```csv
SampleNumber,Barcode1,Barcode2,PCR_product,F_primer,R_primer
Sample_01,AAGTCCAA,GGAGTACT,target_01,ACGT...,TGCA...
Sample_02,GACCTGAA,GGACTTGG,target_02,ACGT...,TGCA...
```

### Validation Rules

- Sample names must not be empty or duplicate.
- Sample names transformed into output filenames must not collide; alphanumeric characters, underscores, and hyphens are recommended.
- `unassigned` is a reserved sample name (case-insensitive).
- `Barcode1` must not be empty; sequences permit only `A/C/G/T/N` (case-insensitive).
- Each barcode must contain at least one non-`N` informational base.
- Single barcode sequences or dual barcode combinations must not be duplicate.
- If any row specifies `Barcode2`, all rows must specify `Barcode2`; mixing single and dual barcodes within the same table is not supported.
- The `Barcode2` column may be omitted entirely or left blank across all rows; both enter single barcode mode.
- Does not automatically reverse-complement `Barcode2`. Enter the observed sequence on the corresponding read.

See the template at [`examples/sample_barcodes.csv`](examples/sample_barcodes.csv). Ultraplex CSV files cannot be directly provided; the repository provides a [conversion script](scripts/seqmux_to_ultraplex_csv.py) to convert SeqMux tables into Ultraplex format for comparison.

---

## Matching and Orientations

### Mode Selection

| Mode | Condition | Target Position |
|:---|:---|:---|
| Dual barcode PE | Table has `Barcode2`, and `-I` is supplied | 5′ ends of both mates |
| Dual barcode SE | Table has `Barcode2`, and `-I` is omitted | 5′ and 3′ ends of the same read |
| Single barcode | Table has no valid `Barcode2` | 5′ end of R1 or SE read |

In PE single barcode mode, R2 is output alongside R1 based on R1 assignment; R2 does not participate in barcode recognition.

### PE Orientations and Canonicalization

| Orientation | Matching Condition |
|:---|:---|
| `canonical` | R1 5′ matches Barcode1, R2 5′ matches Barcode2 |
| `swapped` | R1 5′ matches Barcode2, R2 5′ matches Barcode1 |
| `both` (default) | Simultaneously evaluates both orientations |

Under default `both`, the best valid orientation is selected per sample before comparing across samples. Equal scores from distinct samples are classified as `ambiguous`; ties between orientations of the same sample prioritize `canonical`.

By default, swapped matches exchange the entire R1/R2 records. This does not reverse-complement sequences and does not rewrite original mate numbering in read headers. Use `--no-canonicalize` to retain the input mate order.

Orientation parameters only affect dual-barcode PE; SE dual-barcode does not swap ends or search reverse complements.

### Mismatch Thresholds

Barcode distance is calculated as Hamming mismatch on non-`N` positions; observed `N` bases against expected `A/C/G/T` count as mismatches. Barcode `N` positions serve as UMIs and do not contribute to mismatch scoring.

In dual barcode mode, candidates must satisfy:
- `d1 <= mismatches_1`
- `d2 <= mismatches_2`

Candidate ranking uses `d1 + d2`. SE single barcode mode uses `d1`. These thresholds are per-barcode limits, not total pair mismatch ceilings; barcode indels are not allowed.

Fixed 8 bp, non-UMI, pure `A/C/G/T` dual-barcode PE with `mismatches_1=0`, `mismatches_2=0`, and `min_mismatch_delta=0` uses a dedicated exact fast path. Other configurations use the generic matching path.

### Global Confidence Margin

`--min-mismatch-delta D` requires the best eligible sample to have a distance margin of at least `D` against the closest distinct sample:

1. Identify the best sample satisfying per-end mismatch thresholds.
2. If distinct samples tie for the best eligible match, mark as ambiguous directly.
3. When `D > 0`, search the entire comparable sample pool for the closest distinct sample (which does **not** need to satisfy mismatch thresholds).
4. For dual-barcode PE, take the minimal distance across allowed orientations for each competitor; the candidate does not compete against its own alternative orientation.
5. If the competitor distance is not greater than the best distance, or the distance difference is less than `D`, mark as ambiguous.

For example, if the best eligible sample distance is 1 and the closest competitor distance is 3, `D=2` satisfies the margin (`3 - 1 >= 2`). If the competitor distance is 2, the difference is `1 < 2`, resulting in an ambiguous classification. When no other comparable samples exist, no runner-up is fabricated.

`D=0` preserves legacy tie-breaking semantics. Increasing `D` typically increases rejections but does not guarantee results will change (depending on barcode separation geometry and sequencing error rates). `D` is a distance difference, not error probability or statistical confidence.

---

## Processing Pipeline Order

1. Identify barcodes on raw reads, recording pre-filter matching metrics.
2. For swapped matches, exchange mates according to configuration.
3. Extract barcode UMIs and strip matched barcodes unless `--keep-barcodes`.
4. If a TSO pattern is specified, process the 5′ TSO region of the current SE read or current R1.
5. Apply Phred quality trimming.
6. Apply 3′ adapter trimming for current R1/R2.
7. Enforce minimum length; for PE, if either mate is too short, the entire pair is filtered out.
8. Append barcode UMI to headers, compile post-filter statistics, and write by sample; `--counts-only` skips FASTQ writing.

Unmatched and ambiguous reads do not have barcodes stripped, but still undergo TSO processing, quality trimming, adapter trimming, and length filtering. Consequently, unassigned FASTQs are not guaranteed to be byte-identical copies of input files.

### Adapter Defaults

| Mate | Default 3′ Adapter |
|:---|:---|
| Current R1 / SE | `AGATCGGAAGAGCACACGTCTGAA` |
| Current R2 | `AGATCGGAAGAGCGTCGTG` |

Adapter trimming is enabled by default with a minimum overlap of 3 and a maximum error rate of 0.1. `--no-adapter` disables adapter trimming on both mates; passing an empty string disables trimming on a specific mate (e.g., `--adapter-r2 ''`). When specifying custom adapters, provide sequences matching the respective reads and verify with small test datasets.

The adapter matcher uses seed lookup and dynamic programming paths. Testing covers short overlaps, errors, and indels, but does not claim identical outputs to Cutadapt across all parameter combinations.

### Quality and Length

Quality encoding is parsed as Phred+33. Standard mode supports 5′ and 3′ quality trimming. `--nextseq` uses a NextSeq-style 3′ trimming path (does not apply `--quality-cutoff-5`).

`--min-length 0` does not filter zero-length records. If downstream pipelines reject empty reads, explicitly specify the required minimum length.

---

## UMI and TSO

### Barcode UMIs

```csv
SampleNumber,Barcode1,Barcode2
Sample_01,AAGTNNAA,GGAGNNCT
```

An `N` in a barcode sequence indicates extracting that position as a UMI rather than matching a fixed base. The total barcode length includes UMI positions and is stripped by default. `--keep-barcodes` keeps barcode bases on the read without disabling UMI extraction.

Dual-barcode UMIs are concatenated in `Barcode1 + Barcode2` order. PE barcode UMIs are written to headers of both mates using the `rbc:` tag. The current writer replaces spaces in read headers with underscores; this is not a SAM/BAM `RX:Z:` tag.

### TSO Pattern

```bash
seqmux demux -i reads.fastq.gz -b barcodes.csv -o results --tso-pattern NNNNNIII
```

The pattern accepts only `N` and `I`: `N` is extracted into the UMI, while `I` is ignored and trimmed. TSO processing occurs at the read 5′ end after barcode handling; PE only processes the current R1. Reads shorter than the pattern skip TSO processing without raising an error.

Barcode UMI and TSO UMI extraction occur in separate pipeline steps. If combining them, or if input headers already contain `rbc:` or downstream tools depend on mate ID formats, verify downstream compatibility beforehand.

---

## Common Usage Examples

### Single-end single barcode

```csv
SampleNumber,Barcode1
Sample_01,AAGTCCAA
Sample_02,GACCTGAA
```

```bash
seqmux demux -i reads.fastq.gz -b barcodes_single.csv -o results_se -t 4 --compression-level 1
```

### Single-end dual barcode

```bash
seqmux demux -i reads.fastq.gz -b barcodes_dual.csv -o results_dual_se --no-adapter
```

Barcode2 must be located at the actual 3′ suffix of the read; if additional sequence follows Barcode2, adjust inputs accordingly. Ensure the two barcode regions do not overlap; current implementations do not enforce this constraint during matching.

### Tolerating substitution errors with confidence margin

```bash
seqmux demux \
  -i reads_R1.fastq.gz -I reads_R2.fastq.gz \
  -b barcodes.csv -o results_margin \
  --mismatches-1 1 --mismatches-2 1 \
  --min-mismatch-delta 2 \
  -t 8 --compression-level 1
```

Evaluate these parameters against your barcode Hamming separation and error profile; there is no universally optimal D for all sample sheets.

### QC on the first 200,000 pairs

```bash
seqmux demux \
  -i reads_R1.fastq.gz -I reads_R2.fastq.gz \
  -b barcodes.csv -o qc -p first200k \
  --counts-only --max-reads 200000 -t 8
```

`--counts-only` still executes matching, trimming, and filtering, and outputs the summary TSV; it does not represent a pure barcode matching micro-benchmark that skips all processing. To isolate matching metrics, add `--no-adapter -q 0 --quality-cutoff-5 0 --min-length 0` and inspect pre-filter fields.

### Quality and length filtering

```bash
seqmux demux \
  -i reads_R1.fastq.gz -I reads_R2.fastq.gz \
  -b barcodes.csv -o results_filtered \
  -q 20 --min-length 50 \
  -t 8 --compression-level 1
```

### Fixed orientation and omitting unassigned files

```bash
seqmux demux \
  -i reads_R1.fastq.gz -I reads_R2.fastq.gz \
  -b barcodes.csv -o results_canonical \
  --orientation canonical --discard-unassigned \
  -t 8 --compression-level 1
```

`--discard-unassigned` only suppresses writing unassigned FASTQ files; summary statistics still track them.

---

## Full Parameter Reference

### Input, Output, and Resources

| Parameter | Default | Description |
|:---|:---:|:---|
| `-i, --input` | *Required* | Path to R1 or SE FASTQ |
| `-I, --input2` | *None* | Path to R2 FASTQ; activates PE mode |
| `-b, --barcodes` | *Required* | Sample CSV table with header |
| `-o, --out-dir` | `.` | Output directory, created if needed |
| `-p, --prefix` | `seqmux` | Output file prefix (no path separators allowed) |
| `-t, --threads` | `4` | Worker thread count, must be >= 1; does not equal total process threads |
| `--chunk-reads` | `4096` | Reads/pairs per chunk, must be >= 1 |
| `--compression-level` | `6` | Gzip compression level (`1`–`9`; `1` recommended for fast I/O) |
| `--no-gzip` | *Off* | Write uncompressed `.fastq` files |
| `--counts-only` | *Off* | Skip FASTQ writing; still generates summary TSV |
| `--discard-unassigned` | *Off* | Skip writing unassigned FASTQ files; counts are preserved |
| `--summary` | `<prefix>.summary.tsv` | Custom summary path; parent directory must exist |
| `--force` | *Off* | Overwrite existing output files (does not allow input/output path collisions) |
| `--max-reads` | `0` | Upper limit of reads/pairs to process (`0` for unlimited) |
| `--skip-reads` | `0` | Skip this number of reads/pairs before applying `--max-reads` |
| `--quiet` | *Off* | Suppress progress output and terminal summary; summary file is still written |
| `--log-level` | `info` | Logging verbosity: `error`, `warn`, `info`, `debug`, `trace` |

### Matching and Structure

| Parameter | Default | Description |
|:---|:---:|:---|
| `--mismatches-1` | `0` | Maximum mismatches for Barcode1 informational positions |
| `--mismatches-2` | `0` | Maximum mismatches for Barcode2 informational positions |
| `--min-mismatch-delta` | `0` | Minimum confidence distance margin to closest competitor sample |
| `--orientation` | `both` | PE dual barcode orientation: `both`, `canonical`, `swapped` |
| `--no-canonicalize` | *Off* | Do not swap mates on swapped matches |
| `--keep-barcodes` | *Off* | Retain barcode sequences on reads |
| `--tso-pattern` | *None* | 5′ TSO pattern composed of `N` and `I` |

### Trimming and Filtering

| Parameter | Default | Description |
|:---|:---:|:---|
| `-a, --adapter-r1` | Standard R1 adapter | 3′ adapter sequence for current R1/SE |
| `--adapter-r2` | Standard R2 adapter | 3′ adapter sequence for current R2 |
| `--no-adapter` | *Off* | Disable adapter trimming on both mates |
| `--min-adapter-overlap` | `3` | Minimum adapter overlap, must be >= 1 |
| `--adapter-error-rate` | `0.1` | Maximum error rate, finite number in `0.0..=1.0` |
| `-q, --quality-cutoff-3` | `0` | 3′ Phred+33 quality cutoff |
| `--quality-cutoff-5` | `0` | 5′ quality cutoff in standard mode |
| `--nextseq` | *Off* | Use NextSeq 3′ trimming path |
| `-l, --min-length` | `0` | Minimum post-trimming length; PE drops entire pair if either mate fails |

Run `seqmux --help`, `seqmux demux --help`, or `seqmux validate --help` for full built-in help.

---

## Outputs and Statistics

### File Naming

| Mode | Generated Files |
|:---|:---|
| SE sample | `<prefix>_<sanitized_sample>.fastq.gz` |
| SE unassigned | `<prefix>_unassigned.fastq.gz` |
| PE sample | `<prefix>_<sanitized_sample>_R1.fastq.gz`, `..._R2.fastq.gz` |
| PE unassigned | `<prefix>_unassigned_R1.fastq.gz`, `..._R2.fastq.gz` |
| Summary | `<prefix>.summary.tsv` |

`--no-gzip` changes extensions to `.fastq`. Sample files are created lazily upon actual write; samples with zero matching records may not produce output files. The summary TSV lists only samples with counts > 0, sorted by output label; missing samples in the TSV correspond to count 0.

Each output sample file preserves the relative input order of records/pairs, including normalized swapped pairs. To verify output consistency, compare decompressed FASTQ contents rather than requiring identical gzip file byte hashes.

### Metric Definitions and Stages

The counting unit for SE is reads; for PE it is read pairs (e.g., `total_reads=2,000,000` corresponds to 2,000,000 pairs, not 4,000,000 individual mates).

| Metric | Definition |
|:---|:---|
| `total_reads` | Total reads/pairs processed (excluding skipped records) |
| `assigned` | Reads/pairs passing length filtering and uniquely assigned to a sample |
| `unassigned` | Reads/pairs passing length filtering but unmatched or ambiguous |
| `ambiguous` | Ambiguous reads/pairs before length filtering (legacy semantic) |
| `quality_trimmed` | Reads/pairs trimmed by quality; for PE, either mate trimmed counts once |
| `adapter_trimmed` | Reads/pairs trimmed by adapter; for PE, either mate trimmed counts once |
| `too_short` | Reads/pairs removed by minimum length filtering |
| `five_prime_matched_three_prime_missing` | Reserved compatibility field; not incremented by current workers |
| `orientation_canonical` | Canonical matches at matching stage before length filtering |
| `orientation_swapped` | Swapped matches at matching stage before length filtering |
| `match_rate` | Legacy assignment rate (`assigned / total_reads`, includes length filtering loss) |
| `sample:<label>` | Count of reads/pairs assigned to this sample passing length filtering |
| `matched_before_filter` | Reads/pairs uniquely matched to a sample before length filtering |
| `no_match_before_filter` | Reads/pairs with no match before length filtering (excluding ambiguous) |
| `ambiguous_after_filter` | Ambiguous reads/pairs passing length filtering (subset of `unassigned`) |
| `assignment_rate_before_filter` | Matching-stage assignment fraction (`matched_before_filter / total_reads`) |

Quality trimming, adapter trimming, and length filtering can apply to the same record; these counts are not mutually exclusive categories. PE single-barcode runs may also accumulate canonical counts; do not use this as evidence of dual-barcode orientation validation.

### Mathematical Invariants

The following identities are guaranteed to hold:

```text
total_reads = matched_before_filter + no_match_before_filter + ambiguous
total_reads = assigned + unassigned + too_short
sum(sample counts) = assigned
ambiguous_after_filter <= unassigned
```

For dual-barcode PE, additionally:

```text
orientation_canonical + orientation_swapped = matched_before_filter
```

When `too_short = 0`, matching counts before and after filtering are identical. When length filtering is enabled, legacy `ambiguous` can be greater than `unassigned`. To inspect ambiguous reads within unassigned output, use `ambiguous_after_filter`.

### Example Breakdown

For 8 pairs containing 4 unique matches, 2 ambiguous, and 2 unmatched, where length filtering removes half of each category:

```text
total_reads                    8
assigned                       2
unassigned                     2
ambiguous                      2
too_short                      4
matched_before_filter          4
no_match_before_filter         2
ambiguous_after_filter         1
match_rate                     0.250000
assignment_rate_before_filter  0.500000
```

*(Illustrative example; actual summary TSV format is `metric<TAB>count`.)*

v0.4.1 preserves legacy field names, order, and sample row placements as an exact prefix, appending the four new metrics at the end. Downstream parsers are encouraged to parse metrics by name rather than relying on fixed line counts.

---

## Performance and Resource Usage

### Architecture

```text
Input Reader(s) ──> Bounded Chunk Queue ──> Parallel Workers ──> Ordered Single Writer ──> Output FASTQ.gz
```

For multi-threaded PE input, R1 and R2 are read and decompressed by independent reader threads and assembled into pairs. The writer sequentially handles per-sample output buffering and gzip encoding. `--threads` controls worker thread count, not total process threads.

Task and result channels use bounded capacities. The ordered writer's pending queue for out-of-order chunks does not currently enforce a separate hard memory cap. Tested large amplicon workloads demonstrate stable RSS, but constant memory is not strictly guaranteed under adversarial input conditions.

### Historical Allocator A/B Benchmarks

The [M17A Report](benchmarks/m17a_ab_benchmark/M17A_REPORT.md) records benchmarks on a Linux HPC cluster with 2M pairs, 16-CPU allocation, `-t 12 --compression-level 1`, 1 warm-up run and 5 interleaved iterations:

| Build Configuration | Median Wall Time | Throughput (pairs/s) | Peak RSS |
|:---|---:|---:|---:|
| thin LTO, codegen-units=16, System Allocator | 6.01 s | 332,779 | 125.7 MB |
| thin LTO, codegen-units=16, `mimalloc` | 3.54 s | 564,972 | 467.1 MB |
| fat LTO, codegen-units=1, System Allocator | 6.14 s | 325,733 | 146.1 MB |
| fat LTO, codegen-units=1, `mimalloc` | 3.60 s | 555,556 | 457.6 MB |

Fat LTO provided no throughput advantage in these benchmarks; the release profile retains `lto = "thin"`, `codegen-units = 16`.

On a full-scale dataset of 140,771,720 pairs, historical mimalloc execution recorded 195.23 s (~721,056 pairs/s, 466.1 MB RSS). The comparison baseline of 406.80 s originates from historical v0.3.0 benchmarks and should not be cited as concurrent v0.4.0/v0.4.1 A/B runs. Full-scale sample count parity does not automatically imply byte-identical raw compressed gzip files.

### v0.4.1 Default Parameter Regression

The [v0.4.1 Validation Report](docs/V041_VALIDATION.md) records 2M pairs on an HPC node with `-t 8`, default `compression-level=6`, adapter trimming enabled, 1 warm-up run and 3 interleaved iterations:

| Allocator | v0.4.0 Median Time | v0.4.1 Median Time | v0.4.1 Throughput |
|:---|---:|---:|---:|
| System Allocator | 29.76 s | 29.75 s | 67,227 pairs/s |
| `mimalloc` | 28.95 s | 28.54 s | 70,077 pairs/s |

No throughput regression exceeding 5% was observed. Minor variations should not be construed as confirmed speedups. When compression level, thread count, adapter trimming, disk storage, nodes, or build profiles change, absolute numbers are not directly comparable.

### Practical Guidelines

- When runtime is the primary concern, explicitly set `--compression-level 1` and evaluate output file sizes.
- Use the default system allocator when lower peak RSS is required; benchmark `mimalloc` if memory headroom exists.
- Benchmark `-t 4`, `-t 8`, and `-t 12` on your specific CPU allocation rather than assuming more worker threads are always faster.
- Record inputs, sample sheets, source/binary hashes, full command invocations, toolchain versions, and CPU allocations, using warm-ups and interleaved runs.
- Measure realistic workloads that write FASTQ output; `--counts-only` numbers do not reflect end-to-end FASTQ writing throughput.

### Comparison with Peer Tools

| Tool | Conditions to Control When Comparing with SeqMux |
|:---|:---|
| Ultraplex | CSV/orientation conversions, default quality cutoffs, adapter settings, unassigned file outputs, and dictionary precomputation overhead |
| fqtk | Fixed read structures vs mixed orientations, external trimming requirements, output buffering policies |
| Cutadapt | Matching and trimming rules, multi-step pipeline steps, compression backends, thread allocations, and output directory layouts |

Historical benchmarks in [`docs/ULTRAPLEX_BENCH.md`](docs/ULTRAPLEX_BENCH.md) were conducted across different nodes and workloads and should not be used as an absolute speed ranking under uniform current configurations. Verify per-sample assignment and output semantics before comparing throughput.

---

## Verification and Development

```bash
# Code formatting check
cargo fmt --check

# Clippy linter check
cargo clippy --locked --all-targets --all-features -- -D warnings

# Run all unit, integration, and property tests (standard build)
cargo test --locked --all

# Run all tests with mimalloc feature
cargo test --locked --all --features mimalloc-allocator

# Release builds
cargo build --release --locked
cargo build --release --locked --features mimalloc-allocator

# MSRV verification (requires Rust 1.85.0 toolchain)
cargo +1.85.0 check --locked --all-targets --all-features
```

Tests cover SE/PE modes, paired-end synchronization, output path collisions, barcode margins, adversarial and differential adapter trimming, property tests, parameter boundary enforcement, and statistical invariants.

Large production FASTQ datasets are not tracked in git; the repository includes lightweight test fixtures. v0.4.1 validation records report 100% summary compatibility and byte-for-byte decompressed FASTQ SHA-256 match across 2M pairs for both benchmark libraries; these are historical experiment records, not private datasets executed on every CI run.

Commit `1085fd3` passes CI across Linux and Windows runners with both allocators and verified MSRV 1.85 builds. CI passes do not replace biological validation of sample assignment accuracy on specific user libraries.

---

## Known Limitations

1. **Non-Transactional Output Commits**: Failed runs may leave partial FASTQ files in the output directory. Always use clean output directories, verify exit codes, and check final summary files. Do not assume success merely because output files exist.
2. **Pending Queue Memory Under Severe Imbalance**: In situations with severe chunk processing skew, long read lengths, or custom parameters, pending chunk buffering may increase RSS. Constant memory is not mathematically bounded across arbitrary adversarial inputs.
3. **Path Protection Relies on Normalization**: Path checking primarily compares normalized paths. Scenarios involving hard links or symbolic link aliases require stronger filesystem identity checks. Do not specify input files or their aliases as output paths.
4. **UMI Header Formatting Downstream Compatibility**: Appending `rbc:`, replacing header whitespace with underscores, and preserving original mate numbers alters read headers. Some downstream tools or subsequent `seqmux validate` checks may be sensitive to these formats.
5. **Overlapping SE Dual Barcodes Not Enforced**: For short SE reads, callers must ensure the 5′ and 3′ barcodes do not overlap. This edge case is not validated during matching.
6. **Adapter Trimming Order Relative to Orientation**: Mate canonicalization occurs before adapter trimming. When configuring custom adapters for each mate, verify that the adapter sequences match the canonicalized mate roles.
7. **No Phred+64 Quality Option**: Only Phred+33 quality scores are supported. Do not run quality trimming directly on inputs with unknown or legacy Phred+64 encodings.
8. **File Handle Scaling with Sample Counts**: Sample writer handles are retained until process termination once lazily created. Extremely large sample sheets may encounter OS open file descriptor limits or gzip buffer memory pressures.

---

## Frequently Asked Questions (FAQ)

### Why is `assigned + unassigned` less than `total_reads`?
The difference is `too_short` (reads filtered out by `--min-length`). Inspect `matched_before_filter` to distinguish barcode matching failures from retention loss caused by length filtering.

### Why is `ambiguous` greater than `unassigned`?
`ambiguous` is counted before length filtering, whereas `unassigned` is counted after length filtering. Use `ambiguous_after_filter` to inspect the ambiguous subset remaining within unassigned reads.

### Why did results not change after increasing `min-mismatch-delta`?
The best eligible samples may already be separated from all competing samples by a distance exceeding D, or the comparable sample pool may be small. Check barcode Hamming separation and error distributions rather than assuming the parameter did not take effect.

### Why is the assignment rate lower than expected?
Verify that barcodes are located at the start of reads, that `Barcode2` in the sample sheet matches the observed sequence orientation, that the library is an inline barcode library, and whether orientation is restricted to `canonical`. Compare pre-filter metrics between exact and mismatch-tolerant configurations. Avoid increasing mismatches solely to boost assignment counts.

### Why do output R1 headers still contain `/2` or Illumina mate 2 identifiers?
Default canonicalization swaps entire mate records but preserves read header strings. Output R1 designates the Barcode1-containing mate, which may not correspond to the sequencer's original R1 label. If downstream pipelines require original mate designations, specify `--no-canonicalize` and verify compatibility.

### Why is a specific sample output file missing?
Output files are created lazily upon the first written record. If a sample has zero matches or all matches are eliminated by length filtering, no file is generated, and no non-zero `sample:` line appears in the summary TSV.

### Why does the v0.4.1 summary TSV have four additional rows?
These are appended pre-filter and post-filter metrics added for clarity while preserving existing field order and prefixes. Scripts parsing fixed row counts should be updated to parse by metric name.

### Why is default execution slower than benchmark tables?
The default gzip compression level is 6, whereas many high-throughput benchmarks explicitly use `--compression-level 1`. Differences in allocator, thread allocation, adapter trimming, hardware, and storage systems also affect performance. Compare runs under identical command flags.

### What does `seqmux validate` do?
It validates FASTQ record parsing, sequence and quality string length consistency, and PE record count and normalized read ID synchronization. It outputs `records`, `bases`, `min_length`, `max_length`, `paired`, and `status`. It does not validate barcode tables or calculate demultiplexing accuracy.

### When should `--force` be used?
By default, SeqMux aborts if target output files already exist. Specify `--force` only when intentionally overwriting existing outputs. Using dedicated output directories for each run is recommended to prevent mixing partial results from interrupted runs.

---

## License and Acknowledgements

### License

MIT License — see [LICENSE](LICENSE).

### Acknowledgements

SeqMux builds upon and takes inspiration from outstanding tools and algorithms developed by the NGS bioinformatics and open-source communities:

- **[Ultraplex](https://github.com/ulelab/ultraplex)**: Provided architectural inspiration for unified inline barcode demultiplexing and 3′ adapter trimming workflows.
- **[fqtk](https://github.com/fulcrumgenomics/fqtk)**: Informed high-throughput Rust-based FASTQ stream processing design and confidence margin concepts.
- **[Cutadapt](https://github.com/marcelm/cutadapt)**: Established foundational principles and gold standards for sequencing adapter search and quality trimming algorithms.
- **[needletail](https://github.com/onecodex/needletail)**: High-speed FASTQ parsing backend.
- **[flate2](https://github.com/rust-lang/flate2)** & **[zlib-rs](https://github.com/trifectatechfoundation/zlib-rs)**: Pure-Rust accelerated gzip decompression and compression streaming backends.
- **[mimalloc](https://github.com/microsoft/mimalloc)**: High-performance optional memory allocator.

When reporting issues, please include the SeqMux version (`seqmux --version`), operating system, allocator build, full command invocation, minimal anonymized FASTQ/CSV fixtures, summary TSV, and error logs to facilitate reproduction.
