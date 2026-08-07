# Changelog

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
