# SeqMux vs Ultraplex benchmark

Date: 2026-08-14  
Ultraplex: conda env `PCR` (`ulelab/ultraplex`, Python 3.9)  
SeqMux: `target/release/seqmux` v0.2.0  
Host for 200k smoke: login node; 2M: SLURM `qcpu_23i`

Ultraplex is the historical inspiration for SeqMux, **not** a drop-in replacement.
This document records a fair comparison on the lab I464 dual-barcode PE library.

## Semantic mapping

| Topic | SeqMux | Ultraplex |
|-------|--------|-----------|
| Sample table | `SampleNumber,Barcode1,Barcode2` CSV | Ultraplex CSV: `5p,3p:name` rows |
| Dual PE barcodes | Barcode1 @ R1 **5′**, Barcode2 @ R2 **5′** | 5′ on R1 + 3′ on **RC(R2)** (= R2 5′) |
| Convert lab table | — | `scripts/seqmux_to_ultraplex_csv.py` writes `Barcode1,RC(Barcode2):name` |
| Both orientations | `--orientation both` | emit swapped rows with `--include-swapped`, then collapse `*_swapped` |
| Default mismatches | 0 / 0 | CLI help says 0 / 0 in this install (README historically differed) |
| Default quality trim | off (`-q 0`) | **on** (`-q 30`); set `-q 0` for assignment parity |
| Default adapter trim | **on** (Illumina; `--no-adapter` to disable) | **on** (Illumina); leave default or set `-q 0` and accept adapter work |
| Reference dict | HashMap / Hamming | builds 5^L dict unless `--dont_build_reference` |
| `--ignore_no_match` | `--discard-unassigned` | **avoid** with this build: `merge_defaultdicts` TypeError when values are `None` |

## Result parity (I464 first 200,000 pairs)

Exact mismatch=0, `-q 0`, Ultraplex `--dont_build_reference`.

| Mode | SeqMux assigned | Ultraplex assigned | Per-sample |
|------|----------------:|-------------------:|------------|
| Canonical only | 62,944 (31.47%) | 62,944 | **exact match** |
| Both orientations | 124,867 (62.43%) | 124,867 (collapse `*_swapped`) | **exact match** |

Compare:

```bash
python3 scripts/compare_ultraplex_counts.py \
  tmp/upx_bench/seqmux_can_t8/can.summary.tsv \
  tmp/upx_bench/ultraplex_rc_t8 --prefix upx
```

## Speed (200k pairs, gzip write, login node)

Throughput = pairs / wall seconds.

| Tool | Threads | Wall (s) | Peak RSS | pairs/s | Notes |
|------|--------:|---------:|---------:|--------:|-------|
| SeqMux | 1 | 0.96 | 50 MB | ~208 k | `--discard-unassigned`, gzip lvl 1, no adapter |
| SeqMux | 8 | 0.68 | 64 MB | ~294 k | same |
| SeqMux | 8 | 0.68 | 75 MB | ~294 k | also write unassigned |
| Ultraplex | 1 | 34.69 | 58 MB | ~5.8 k | default adapters, `-q 0`, `--dont_build_reference` |
| Ultraplex | 8 | 6.89 | 53 MB | ~29 k | same, canonical CSV |
| Ultraplex | 8 | 8.21 | 58 MB | ~24 k | both-orientation CSV (70 rows) |

Approx wall-clock speedup (SeqMux / Ultraplex): **~36×** at 1 thread, **~10×** at 8 threads on this 200k write workload.

Without `--dont_build_reference`, Ultraplex spends minutes building 5^8 dictionaries per worker for 8-mer barcodes (not usable for this library).

## Scale check (2,000,000 pairs, SLURM `qcpu_23i`, 16 threads)

Job `2314929`. Canonical write path; SeqMux `--discard-unassigned` + gzip lvl 1; Ultraplex `-q 0 --dont_build_reference` (still writes `*_no_match*` files).

| Tool | Wall (s) | Peak RSS | assigned | Per-sample |
|------|---------:|---------:|---------:|------------|
| SeqMux | 4.55 | 76 MB | 629,825 | — |
| Ultraplex | 21.54 | 59 MB | 629,825 | **exact match** |
| SeqMux both `--counts-only` | 4.96 | 17 MB | (reference) | — |

Speedup at 2M / 16 threads ≈ **4.7×** (Ultraplex still does adapter trim + writes all no-match files).

## Caveats

1. Ultraplex still runs Cutadapt-style adapter trimming by default; SeqMux does not unless `-a` is set. That adds CPU to Ultraplex beyond demux.
2. Ultraplex always writes gzip sample FASTQs (plus many `*_no_match*` files unless `--ignore_no_match`, which crashes here). SeqMux can `--counts-only`.
3. Ultraplex CSV cannot express SeqMux “both orientations” natively; the swapped-row workaround doubles the barcode table.
4. Lab production demux remains the Python exact matcher (parity already proven in `docs/COMPATIBILITY.md`). Ultraplex is a third-party speed/result reference.

## Reproduce

```bash
# convert + extract
python3 scripts/seqmux_to_ultraplex_csv.py \
  tests/fixtures/I464-469erdai_barcode_and_name.csv \
  -o tmp/upx_bench/i464.ultraplex.csv
python3 scripts/extract_pe_fastq_head.py \
  -i /path/to/I464_1.fq.gz -I /path/to/I464_2.fq.gz \
  -o tmp/upx_bench/i464_200k_1.fq.gz -O tmp/upx_bench/i464_200k_2.fq.gz -n 200000

# SeqMux
./target/release/seqmux demux \
  -i tmp/upx_bench/i464_200k_1.fq.gz -I tmp/upx_bench/i464_200k_2.fq.gz \
  -b tests/fixtures/I464-469erdai_barcode_and_name.csv \
  -o tmp/upx_bench/seqmux_can -p can --orientation canonical -t 8 \
  --discard-unassigned --compression-level 1 --force

# Ultraplex (PCR env)
ultraplex -i tmp/upx_bench/i464_200k_1.fq.gz -i2 tmp/upx_bench/i464_200k_2.fq.gz \
  -b tmp/upx_bench/i464.ultraplex.csv -d tmp/upx_bench/ultraplex_rc/ \
  -o upx -t 8 -m5 0 -m3 0 -q 0 -l 1 -ig --dont_build_reference
```

SLURM 2M-pair job: `tmp/upx_bench/run_ultraplex_bench.slurm`.
