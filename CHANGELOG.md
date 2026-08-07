# Changelog

## 0.1.0 — 2026-08-07

Initial public release of **SeqMux**.

- Single-end and paired-end FASTQ demultiplex
- 5′ / linked 3′ barcodes, UMI (`N`) extraction to `rbc:` header tag
- Hamming mismatch matching with ambiguous-tie → unassigned
- Quality trim (BWA-style) and 3′ adapter trim
- `--three-prime-only` and TSO pattern support
- Chunked multi-threaded pipeline with ordered writer
- `seqmux validate` for FASTQ / pair checks
- Linux + Windows CI
