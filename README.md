# SeqMux

**SeqMux** is a fast, portable FASTQ demultiplexer written in Rust.

It is inspired by [Ultraplex](https://github.com/ulelab/ultraplex), redesigned as a **single binary** with streaming I/O, no Python/Conda/pigz/SLURM dependency, and identical behaviour on Linux and Windows.

## Features (v0.1)

- Single-end and paired-end FASTQ (`.fastq` / `.fq` / `.gz`)
- 5′ barcode demultiplex with optional linked 3′ barcodes
- `N` positions as UMI (moved into the read header as `rbc:`)
- Hamming-distance barcode matching with mismatch tolerance and tie → unassigned
- Quality trimming (BWA-style) and 3′ adapter trimming
- `--three-prime-only` paired-end mode
- TSO pattern (`N` = UMI, `I` = trim-only)
- Multi-threaded chunked pipeline with **ordered** per-sample output
- Gzip output without shelling out to pigz
- Summary TSV + stderr report

## Install

### From source

```bash
git clone https://github.com/Caizhaohui/SeqMux.git
cd SeqMux
cargo build --release
# binary: target/release/seqmux
```

### Requirements

- Rust stable (1.74+)
- No external runtime tools

## Quick start

**Single-end:**

```bash
seqmux demux \
  -i reads.fastq.gz \
  -b barcodes.csv \
  -o results \
  -t 8
```

**Paired-end:**

```bash
seqmux demux \
  -i reads_R1.fastq.gz \
  -I reads_R2.fastq.gz \
  -b barcodes.csv \
  -o results \
  -t 8
```

**Validate FASTQ:**

```bash
seqmux validate -i reads.fastq.gz
seqmux validate -i r1.fastq.gz -I r2.fastq.gz
```

## Barcode CSV

Ultraplex-compatible format:

```csv
NNNATGNN:sample1,
NNNCCGNN,ATG:sample2,TCA:sample3
NNNCACNN,
```

Rules:

- First column = 5′ barcode (or 3′-only barcodes when using `--three-prime-only`)
- Additional columns = linked 3′ barcodes
- Optional sample names after `:`
- `N` = UMI base (not used for identity matching)
- Only `A/C/G/T/N` allowed

## Common options

| Option | Description |
|--------|-------------|
| `-i / --input` | R1 or single-end FASTQ |
| `-I / --input2` | R2 FASTQ |
| `-b / --barcodes` | Barcode CSV |
| `-o / --out-dir` | Output directory |
| `-p / --prefix` | Output prefix (default `seqmux`) |
| `-t / --threads` | Worker threads |
| `--mismatches-5` | Allowed 5′ mismatches |
| `--mismatches-3` | Allowed 3′ mismatches |
| `--keep-barcodes` | Do not trim barcode bases from sequence |
| `--discard-unassigned` | Drop unassigned reads |
| `-a / --adapter-r1` | 3′ adapter for R1 |
| `-q / --quality-cutoff-3` | 3′ Phred quality cutoff |
| `-l / --min-length` | Minimum length after trim |
| `--force` | Overwrite existing outputs |
| `--no-gzip` | Write plain FASTQ |

Run `seqmux demux --help` for the full list.

## Output layout

Single-end:

```text
<prefix>_<sample>.fastq.gz
<prefix>_unassigned.fastq.gz
<prefix>.summary.tsv
```

Paired-end:

```text
<prefix>_<sample>_R1.fastq.gz
<prefix>_<sample>_R2.fastq.gz
```

Without sample names, files use barcode labels such as `5bc_NNNATGNN` or `5bc_…_3bc_…`.

## Design notes

SeqMux **does not** build Ultraplex-style `5^L` reference dictionaries. Barcodes are compiled to informative positions and matched with Hamming distance (exact hash planned as a later optimisation). Processing uses:

```text
Reader thread → bounded chunks → N workers → ordered writer
```

Workers never open output files; the writer alone finishes gzip streams. There are no temporary `_tmp_thread_*` files and no final `cat` merge.

## Development

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all
cargo build --release
```

See [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) for the full milestone plan (originally written for the Ultraplex → Rust rewrite, applied here under the SeqMux name).

## License

MIT — see [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md).

## Citation / inspiration

If you use SeqMux in published work, please also cite Ultraplex where appropriate:

- Wilkins et al., Ultraplex — https://github.com/ulelab/ultraplex
