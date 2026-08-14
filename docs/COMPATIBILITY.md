# Compatibility with the lab Python demux

Reference implementations:

- `/hpcfs/fhome/caizhh/03_PCR_analysis/25_464469erdaimix/scripts/demux_i464_fastq.py`
- `/hpcfs/fhome/caizhh/03_PCR_analysis/23_I395erdaimix/scripts/demux_i395_fastq.py`

Both register exact 8-mer pairs in **two orientations** and trim 8 bp from each 5′ end. No mismatch, no quality trim.

## Counts (must match)

With `--orientation both` (default), `--mismatches-1 0`, `--mismatches-2 0`:

- total pairs
- assigned / unassigned
- per-sample assigned counts
- orientation split (`orientation_canonical` = Python `R1_Barcode1_R2_Barcode2`; `orientation_swapped` = Python `R1_Barcode2_R2_Barcode1`)

must equal the Python QC tables. Gzip bytes need not match.

## Documented differences (FASTQ contents, not assignment)

| Topic | Python | SeqMux default |
|-------|--------|----------------|
| Output R1 | mate that carried **Barcode2** | mate that carried **Barcode1** (`--no-canonicalize` keeps original R1/R2) |
| Read header `/1` `/2` | follows the swapped physical read | same: header travels with the sequence |
| Unassigned FASTQ | not written | written unless `--discard-unassigned` |
| UMI tag `rbc:` | not used (barcodes have no `N`) | added only when barcode patterns contain `N` |
| Quality / adapter trim | none | off by default; if enabled, runs **after** barcode match |

## Not compared

- Python `unmatched_barcode_counts.csv` top-N histogram (SeqMux does not emit this yet)
- fastp merge / mutation counting downstream of demux
