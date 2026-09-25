# Changelog

## 0.3.0 — 2026-09-25

Release candidate. Not tagged.

### Changed

- 3′ adapter search keeps stride-1 seed completeness and uses a fixed packed-8-mer lookup plus short-overlap candidate filtering. Generic DP remains the correctness reference.
- Paired-end inputs with `-t > 1` decompress R1 and R2 on two reader threads. One ordered writer is unchanged.
- Recommended I395/I464 production settings: **`-t 12`**, **`--compression-level 1`**, mismatch 0, adapter trimming on, `--orientation both` with canonicalization. The CLI compression default stays 6.

### Validated

Full gzip-in / gzip-level-1-out runs on a 16-CPU allocation. Per-sample assignment matched the lab Python exact reference for all 35 I464 samples and all 18 I395 samples. Gzip streams were integrity-checked; they were not compared byte-for-byte with another artifact.

I464: 140,771,720 pairs; wall 406.80 s (6 min 47 s); 346,047 pairs/s; assigned 87,477,095; unassigned 53,294,625; canonical 45,393,500; swapped 42,083,595; adapter trimmed 1,240,406; peak RSS 155,860 kB; output 10,901,284,829 B.

I395: 60,523,000 pairs; wall 177.03 s (2 min 57 s); 341,880 pairs/s; assigned 39,455,648; unassigned 21,067,352; canonical 18,606,371; swapped 20,849,277; adapter trimmed 795,256; peak RSS 128,228 kB; output 5,168,089,121 B.

Checked on representative samples: `gzip -t`, paired record counts, paired read IDs, head/tail record order against the input, canonical Barcode1-on-R1, and SHA-256 of those files.

## 0.2.1 — 2026-08-14

### Changed

- 3′ Illumina adapter trimming is **on by default** (R1 `AGATCGGAAGAGCACACGTCTGAA`, R2 `AGATCGGAAGAGCGTCGTG`, same as Ultraplex defaults). Use `--no-adapter` to disable; empty `-a` / `--adapter-r2` disables that mate only.
- Dual-barcode PE `--orientation both` remains the default (unchanged since 0.1.2).

### Docs

- `docs/ULTRAPLEX_BENCH.md`: SeqMux vs Ultraplex speed/assignment benchmark on I464
- Helper scripts for Ultraplex CSV conversion and count comparison

## 0.2.0 — 2026-08-13

### Added

- `--max-reads N` / `--skip-reads N` for smoke and shard runs
- `--counts-only` for assignment QC without writing FASTQ
- Progress on stderr: processed pairs, k pairs/s, assigned %, elapsed
- Summary `match_rate`
- `docs/REAL_DATA.md`, `docs/COMPATIBILITY.md`, `docs/MISMATCH_QC.md`
- Helper scripts: `scripts/compare_python_counts.py`, `scripts/audit_mismatch1.py`

### Changed

- Dual-barcode PE matches **both insert orientations** by default (from 0.1.2)
- Pipeline order: **barcode match → trim barcodes → quality → adapter**
- Gzip via `flate2` + `zlib-rs` (pure Rust); release profile uses thin LTO
- Validated on full I464 (140.8M) and I395 (60.5M) vs lab Python demux counts

### Docs

- Mismatch=1 stays opt-in (default 0); I464 full QC recorded in `docs/MISMATCH_QC.md`

## 0.1.2 — 2026-08-13

### Fixed

- **Dual-barcode PE now matches both insert orientations by default.**  
  Real amplicon libraries (I395 / I464) have ~50/50 `Barcode1@R1` vs `Barcode2@R1`.  
  v0.1.1 only accepted the canonical orientation and silently dropped half of assignable pairs.
- Swapped-orientation matches are **canonicalized** so output R1 carries Barcode1 (use `--no-canonicalize` to keep original mates).
- Summary reports `orientation_canonical` / `orientation_swapped`.

### CLI

- `--orientation both|canonical|swapped` (default `both`)
- `--no-canonicalize`

## 0.1.1 — 2026-08-07

### Breaking

- **Barcode input is now SeqMux sample table CSV only**  
  Required header: `SampleNumber`, `Barcode1`, optional `Barcode2`.  
  Ultraplex-style CSV (`NNNATGNN:sample,...`) is **no longer supported**.
- Dual-barcode PE matching: `Barcode1` @ R1 5′, `Barcode2` @ R2 5′
- Dual-barcode SE matching: `Barcode1` @ R1 5′, `Barcode2` @ R1 3′
- CLI: `--mismatches-1` / `--mismatches-2` (old `--mismatches-5`/`-3` kept as aliases)
- Removed Ultraplex-specific `--three-prime-only` mode

### Docs

- README rewritten around the sample barcode table format
- Added `examples/sample_barcodes.csv`

## 0.1.0 — 2026-08-07

Initial public release of **SeqMux**.

- Single-end and paired-end FASTQ demultiplex
- Chunked multi-threaded pipeline with ordered writer
- `seqmux validate` for FASTQ / pair checks
- Linux + Windows CI
