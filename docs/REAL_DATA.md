# Real experimental data (not in git)

FASTQ files stay on the lab filesystem. SeqMux only vendors a 5-pair fixture under `tests/fixtures/i464_real_*.fastq`.

## I464-469erdaimix

| Item | Path |
|------|------|
| R1 | `/hpcfs/fhome/caizhh/03_PCR_analysis/25_464469erdaimix/E260702007_L01_464469erdaimix_1.fq.gz` |
| R2 | `/hpcfs/fhome/caizhh/03_PCR_analysis/25_464469erdaimix/E260702007_L01_464469erdaimix_2.fq.gz` |
| Table | `/hpcfs/fhome/caizhh/03_PCR_analysis/25_464469erdaimix/I464-469erdai_barcode_and_name.csv` |
| Python QC | `.../results/barcode_demux_qc.csv`, `sample_read_counts.csv`, `barcode_orientation_qc.csv` |

### Python exact demux (reference)

| Metric | Value |
|--------|-------|
| total pairs | 140,771,720 |
| matched | 87,477,095 |
| unmatched | 53,294,625 |
| match rate | 0.621411 |
| R1=BC1 | 45,393,500 |
| R1=BC2 | 42,083,595 |

### SeqMux v0.2.0 full counts-only（历史，2026-08-13）

| Metric | SeqMux | vs Python |
|--------|--------|-----------|
| total | 140,771,720 | exact |
| assigned | 87,477,095 | exact |
| unassigned | 53,294,625 | exact |
| orientation_canonical | 45,393,500 | exact |
| orientation_swapped | 42,083,595 | exact |
| per-sample (35) | — | **all match** |
| wall clock | 5:58 | — |
| peak RSS | 16 MB | — |
| throughput | ~393 k pairs/s | — |

Compare:

```bash
python3 scripts/compare_python_counts.py \
  tmp/m12_i464/i464.summary.tsv \
  /hpcfs/fhome/caizhh/03_PCR_analysis/25_464469erdaimix/results/sample_read_counts.csv \
  --expect-total 140771720 --expect-assigned 87477095 \
  --expect-canonical 45393500 --expect-swapped 42083595
```

Smoke (does not write FASTQ):

```bash
seqmux demux \
  -i E260702007_L01_464469erdaimix_1.fq.gz \
  -I E260702007_L01_464469erdaimix_2.fq.gz \
  -b I464-469erdai_barcode_and_name.csv \
  -o /tmp/seqmux_i464_smoke \
  --counts-only --max-reads 200000 -t 8 --force
```

## I395erdaimix

| Item | Path |
|------|------|
| R1 | `/hpcfs/fhome/caizhh/03_PCR_analysis/23_I395erdaimix/E260617005_L01_395erdaimix_1.fq.gz` |
| R2 | `/hpcfs/fhome/caizhh/03_PCR_analysis/23_I395erdaimix/E260617005_L01_395erdaimix_2.fq.gz` |
| Table | `/hpcfs/fhome/caizhh/03_PCR_analysis/23_I395erdaimix/I395erdai_gst_barcode_and_name.csv` |

### Python / SeqMux (exact match)

| Metric | Value |
|--------|-------|
| total pairs | 60,523,000 |
| matched | 39,455,648 |
| unmatched | 21,067,352 |
| match rate | 0.651912 |
| R1=BC1 | 18,606,371 |
| R1=BC2 | 20,849,277 |
| SeqMux wall | 2:32 |
| peak RSS | 14 MB |
| per-sample (18) | **all match** |

## SeqMux v0.3 production (2026-09-25)

Same inputs and barcode tables. Command: mismatch 0 (default), adapter trimming on (default), `--compression-level 1`, **`-t 12`**, 16-CPU `qcpu_18i` allocation. Chosen because median wall time on the 2M I464 production path was 5.85 s (`-t 12`) versus 6.09 s (`-t 8`) and 5.87 s (`-t 16`).

| Metric | I464 | I395 |
|--------|------|------|
| total pairs | 140,771,720 | 60,523,000 |
| wall | 406.80 s (6 min 47 s) | 177.03 s (2 min 57 s) |
| throughput | 346,047 pairs/s | 341,880 pairs/s |
| user CPU | 2671.28 s | 1155.31 s |
| sys CPU | 1340.83 s | 590.38 s |
| CPU % | 986% | 986% |
| peak RSS | 155,860 kB | 128,228 kB |
| assigned | 87,477,095 | 39,455,648 |
| unassigned | 53,294,625 | 21,067,352 |
| adapter_trimmed | 1,240,406 | 795,256 |
| output size | 10,901,284,829 B | 5,168,089,121 B |
| per-sample vs Python | all 35 samples | all 18 samples |
| orientation canonical / swapped | 45,393,500 / 42,083,595 | 18,606,371 / 20,849,277 |

CLI `--compression-level` 的默认值仍是 6。上表使用推荐的高吞吐设置 `--compression-level 1`，没有改默认值。

核对过的输出性质：`gzip -t`、R1/R2 记录数、mate ID、头尾记录相对输入的顺序、canonical 方向（Barcode1 在输出 R1），以及代表性文件的 SHA-256。分配计数与 Python 参考逐样本一致。没有把整套 gzip 输出做成与另一份产物逐字节相同的测试。

## Performance notes（v0.2.0 历史，不是当前生产数字）

| Workload | Wall | RSS | Notes |
|----------|------|-----|-------|
| I464 counts-only (140.8M) | 5:58 | 16 MB | gzip decompress bound; `zlib-rs` backend |
| I395 counts-only (60.5M) | 2:32 | 14 MB | same |
| I464 write 2M pairs (gzip lvl 1, discard unassigned) | 5.8 s | 77 MB | 70 sample FASTQ.gz files, ~84 MB output |

Login-node runs may be SIGKILL'd after ~4 min; use SLURM (`qcpu_23i`) for full datasets.

## Notes

- `--orientation canonical` on I464 200k subsample assigns 62,944 / 200,000 (31.47%), matching v0.1.1.
- `--orientation both` (default) assigns 124,867 / 200,000 (62.43%), matching Python.
- SeqMux output R1 = Barcode1 mate after canonicalize; Python scripts wrote Barcode2 as R1. See `COMPATIBILITY.md`.
- Mismatch=1 QC: see `MISMATCH_QC.md`.
