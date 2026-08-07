# SeqMux

**SeqMux** is a fast, portable FASTQ demultiplexer written in Rust.

It demultiplexes single-end and paired-end FASTQ using a **sample barcode table** (dual or single barcode), with quality/adapter trimming and multi-threaded streaming I/O. No Python, Conda, pigz, or SLURM required.

## Features (v0.1)

- Single-end and paired-end FASTQ (`.fastq` / `.fq` / `.gz`)
- **SeqMux sample table CSV** (`SampleNumber`, `Barcode1`, `Barcode2`, …)
- Dual-barcode demux (PE: R1/R2 5′; SE: R1 5′ + R1 3′)
- Single-barcode demux (Barcode1 only at R1 5′)
- Optional `N` bases in barcodes as UMI (written to header as `rbc:`)
- Hamming-distance matching with mismatch tolerance; ties → unassigned
- Quality trimming (BWA-style) and 3′ adapter trimming
- Multi-threaded chunked pipeline with **ordered** per-sample output
- Gzip output without shelling out to external tools
- Summary TSV + stderr report

## Install

```bash
git clone https://github.com/Caizhaohui/SeqMux.git
cd SeqMux
cargo build --release
# binary: target/release/seqmux
```

Requirements: Rust stable (1.74+). No external runtime tools.

## Quick start

**Paired-end (recommended for dual barcode):**

```bash
seqmux demux \
  -i mix_R1.fastq.gz \
  -I mix_R2.fastq.gz \
  -b sample_barcodes.csv \
  -o results \
  -t 8
```

**Single-end:**

```bash
seqmux demux \
  -i reads.fastq.gz \
  -b sample_barcodes.csv \
  -o results \
  -t 8
```

**Validate FASTQ:**

```bash
seqmux validate -i reads.fastq.gz
seqmux validate -i r1.fastq.gz -I r2.fastq.gz
```

## Sample barcode table (required format)

SeqMux reads a **CSV with a header row**. This is the only supported barcode input.

### Minimal example

```csv
SampleNumber,Barcode1,Barcode2
I464469-A1,AAGTCCAA,GGAGTACT
I464469-A2,GACCTGAA,GGACTTGG
I464469-A3,GGCTTAAG,TGGATCGA
```

### Full example (extra columns are ignored)

Same layout as experimental tables such as `I464-469erdai_barcode_and_name.csv`:

```csv
SampleNumber,Barcode1,Barcode2,PCR_product,rawdata1,rawdata2,library_round,F_primer,R_primer
I464469-A1,AAGTCCAA,GGAGTACT,atgcca...,mix_1.fq.gz,mix_2.fq.gz,A,fwd...,rev...
I464469-A2,GACCTGAA,GGACTTGG,agcgcg...,mix_1.fq.gz,mix_2.fq.gz,A,fwd...,rev...
```

A trimmed example file is included at [`examples/sample_barcodes.csv`](examples/sample_barcodes.csv).

### Required columns

| Column | Required | Description |
|--------|----------|-------------|
| `SampleNumber` | **Yes** | Sample ID / output name (unique) |
| `Barcode1` | **Yes** | First barcode (DNA: A/C/G/T/N) |
| `Barcode2` | Optional* | Second barcode for dual-index demux |

\*If any row has `Barcode2`, **all** rows must have `Barcode2`.  
If `Barcode2` is omitted entirely, SeqMux runs in **single-barcode** mode.

### Column name aliases (case-insensitive)

| Logical field | Accepted headers |
|---------------|------------------|
| Sample | `SampleNumber`, `Sample`, `Sample_Name`, `Sample_ID`, `Name` |
| Barcode1 | `Barcode1`, `Barcode_1`, `BC1`, `Index1`, `i7` |
| Barcode2 | `Barcode2`, `Barcode_2`, `BC2`, `Index2`, `i5` |

Extra columns (`PCR_product`, `rawdata1`, `rawdata2`, `library_round`, `F_primer`, `R_primer`, …) are **ignored** and may remain in the file for laboratory records.

### Rules

- Header row is **required**
- Sample names must be unique
- Barcode combinations `(Barcode1, Barcode2)` must be unique
- Bases: only `A/C/G/T/N` (auto-uppercased); `N` = UMI position
- Empty lines are skipped via the CSV reader

### How barcodes are matched

| Mode | When | Matching |
|------|------|----------|
| **Dual barcode + PE** | `Barcode2` present and `-I` given | `Barcode1` @ R1 5′, `Barcode2` @ R2 5′ |
| **Dual barcode + SE** | `Barcode2` present, no `-I` | `Barcode1` @ R1 5′, `Barcode2` @ R1 3′ |
| **Single barcode** | no `Barcode2` column / all empty | `Barcode1` @ R1 5′ only |

After assignment (unless `--keep-barcodes`):

- PE dual: trim Barcode1 from R1 5′ and Barcode2 from R2 5′
- SE dual: trim Barcode1 from 5′ and Barcode2 from 3′
- Single: trim Barcode1 from 5′

Mismatch thresholds: `--mismatches-1` / `--mismatches-2`. Best unique score wins; ties go to `unassigned`.

## Common options

| Option | Description |
|--------|-------------|
| `-i / --input` | R1 or single-end FASTQ |
| `-I / --input2` | R2 FASTQ |
| `-b / --barcodes` | Sample barcode table CSV |
| `-o / --out-dir` | Output directory |
| `-p / --prefix` | Output prefix (default `seqmux`) |
| `-t / --threads` | Worker threads |
| `--mismatches-1` | Allowed mismatches for Barcode1 |
| `--mismatches-2` | Allowed mismatches for Barcode2 |
| `--keep-barcodes` | Do not trim barcode bases |
| `--discard-unassigned` | Drop unassigned reads |
| `-a / --adapter-r1` | 3′ adapter for R1 |
| `--adapter-r2` | 3′ adapter for R2 |
| `-q / --quality-cutoff-3` | 3′ Phred quality cutoff |
| `-l / --min-length` | Minimum length after trim |
| `--force` | Overwrite existing outputs |
| `--no-gzip` | Write plain FASTQ |

Run `seqmux demux --help` for the full list.

## Output layout

Single-end:

```text
<prefix>_<SampleNumber>.fastq.gz
<prefix>_unassigned.fastq.gz
<prefix>.summary.tsv
```

Paired-end:

```text
<prefix>_<SampleNumber>_R1.fastq.gz
<prefix>_<SampleNumber>_R2.fastq.gz
<prefix>_unassigned_R1.fastq.gz
<prefix>_unassigned_R2.fastq.gz
```

Example: sample `I464469-A1` → `seqmux_I464469-A1_R1.fastq.gz`.

## Design notes

- Pipeline: `Reader → bounded chunks → N workers → ordered writer`
- Workers never open output files; writer alone finishes gzip streams
- No temporary `_tmp_thread_*` files and no final `cat` merge

## Development

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all
cargo build --release
```

## License

MIT — see [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md).

## Acknowledgements

Pipeline and trimming ideas were informed by demultiplexing practice in the NGS community, including Ultraplex.
