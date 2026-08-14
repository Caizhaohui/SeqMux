# Changelog

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
