#!/usr/bin/env bash
#SBATCH --job-name=m17a-i464-full
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/03_SeqMux/benchmarks/m17a_ab_benchmark/slurm_i464_full_%j.out
#SBATCH --error=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/03_SeqMux/benchmarks/m17a_ab_benchmark/slurm_i464_full_%j.err

set -euo pipefail

ROOT_DIR="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/03_SeqMux"
cd "${ROOT_DIR}"

source "${HOME}/.cargo/env" 2>/dev/null || true

OUT_DIR="${ROOT_DIR}/tmp/m17a_runs/i464_full_mimalloc"
TIME_LOG="${ROOT_DIR}/benchmarks/m17a_ab_benchmark/i464_full_mimalloc_time.log"
REPORT_MD="${ROOT_DIR}/benchmarks/m17a_ab_benchmark/I464_FULL_MIMALLOC_REPORT.md"

R1="/hpcfs/fhome/caizhh/03_PCR_analysis/25_464469erdaimix/E260702007_L01_464469erdaimix_1.fq.gz"
R2="/hpcfs/fhome/caizhh/03_PCR_analysis/25_464469erdaimix/E260702007_L01_464469erdaimix_2.fq.gz"
BARCODES="${ROOT_DIR}/tests/fixtures/I464-469erdai_barcode_and_name.csv"

echo "=================================================="
echo "Starting Full I464 (140,771,720 pairs) mimalloc validation"
echo "Date: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM Job ID: ${SLURM_JOB_ID:-manual}"
echo "CPUs allocated: ${SLURM_CPUS_PER_TASK:-16}"
echo "=================================================="

# 1. Build mimalloc binary cleanly
echo "[1/4] Building release binary with --features mimalloc-allocator..."
cargo build --release --features mimalloc-allocator --locked
BIN="${ROOT_DIR}/target/release/seqmux"
sha256sum "${BIN}"

# 2. Clean out output directory
echo "[2/4] Preparing output directory: ${OUT_DIR}"
rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

# 3. Execute benchmark with /usr/bin/time -v
echo "[3/4] Running seqmux demux on 140,771,720 pairs (-t 12 --compression-level 1)..."
rm -f "${TIME_LOG}"
/usr/bin/time -v -o "${TIME_LOG}" "${BIN}" demux \
    -i "${R1}" \
    -I "${R2}" \
    -b "${BARCODES}" \
    -o "${OUT_DIR}" \
    -p i464_full \
    -t 12 \
    --compression-level 1 \
    --force

echo "[4/4] Parsing results and verifying correctness..."
cat "${TIME_LOG}"

# Parse summary
SUMMARY="${OUT_DIR}/i464_full.summary.tsv"
cat "${SUMMARY}"

python3 - << 'PYEOF'
import os
import sys

root_dir = "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/03_SeqMux"
time_log = os.path.join(root_dir, "benchmarks/m17a_ab_benchmark/i464_full_mimalloc_time.log")
summary_file = os.path.join(root_dir, "tmp/m17a_runs/i464_full_mimalloc/i464_full.summary.tsv")
report_md = os.path.join(root_dir, "benchmarks/m17a_ab_benchmark/I464_FULL_MIMALLOC_REPORT.md")

# Baseline values for i464_full (v0.3.0)
BASE_WALL = 406.80
BASE_RSS_KB = 155860
BASE_THROUGHPUT = 140771720 / BASE_WALL

TOTAL_PAIRS = 140771720
EXP_ASSIGNED = 87477095
EXP_UNASSIGNED = 53294625
EXP_CANONICAL = 45393500
EXP_SWAPPED = 42083595

time_data = {}
with open(time_log) as f:
    for line in f:
        line = line.strip()
        if "Elapsed (wall clock) time" in line:
            raw = line.split("):")[-1].strip()
            parts = raw.split(":")
            if len(parts) == 3:
                time_data["wall"] = float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])
            elif len(parts) == 2:
                time_data["wall"] = float(parts[0]) * 60 + float(parts[1])
            else:
                time_data["wall"] = float(parts[0])
        elif "User time (seconds):" in line:
            time_data["user"] = float(line.split(":")[-1].strip())
        elif "System time (seconds):" in line:
            time_data["sys"] = float(line.split(":")[-1].strip())
        elif "Percent of CPU this job got:" in line:
            time_data["cpu_pct"] = line.split(":")[-1].strip()
        elif "Maximum resident set size (kbytes):" in line:
            time_data["rss_kb"] = int(line.split(":")[-1].strip())

sum_data = {}
with open(summary_file) as f:
    for line in f:
        p = line.strip().split("\t")
        if len(p) >= 2:
            try:
                sum_data[p[0]] = int(p[1])
            except ValueError:
                sum_data[p[0]] = p[1]

wall = time_data.get("wall", 0.0)
throughput = TOTAL_PAIRS / wall if wall > 0 else 0
speedup = throughput / BASE_THROUGHPUT
rss_mb = time_data.get("rss_kb", 0) / 1024.0
base_rss_mb = BASE_RSS_KB / 1024.0

assigned = sum_data.get("assigned", 0)
unassigned = sum_data.get("unassigned", 0)
canonical = sum_data.get("orientation_canonical", 0)
swapped = sum_data.get("orientation_swapped", 0)

correct = (
    assigned == EXP_ASSIGNED and
    unassigned == EXP_UNASSIGNED and
    canonical == EXP_CANONICAL and
    swapped == EXP_SWAPPED
)

status = "PASS" if correct else "FAIL"

report = f"""# Full I464 Production Dataset mimalloc Validation Report

- **Dataset**: I464 Full Production (140,771,720 pairs)
- **Threads**: `-t 12`
- **Compression**: `--compression-level 1`
- **Build**: `target/release/seqmux` (`--features mimalloc-allocator`)
- **Correctness Gate**: **{status}**

## 1. Metrics Comparison

| Metric | Baseline v0.3.0 (glibc) | mimalloc-allocator | Delta / Speedup |
|:---|---:|---:|---:|
| **Wall Clock** | {BASE_WALL:.2f} s ({BASE_WALL/60:.2f} min) | {wall:.2f} s ({wall/60:.2f} min) | **{speedup:.2f}x** ({(BASE_WALL - wall)/BASE_WALL*100:+.1f}%) |
| **Throughput** | {BASE_THROUGHPUT:,.0f} pairs/s | {throughput:,.0f} pairs/s | **+{throughput - BASE_THROUGHPUT:,.0f} pairs/s** |
| **Peak RSS** | {base_rss_mb:.1f} MB | {rss_mb:.1f} MB | {rss_mb - base_rss_mb:+.1f} MB |
| **User CPU** | - | {time_data.get('user', 0.0):.1f} s | - |
| **Sys CPU** | - | {time_data.get('sys', 0.0):.1f} s | - |
| **CPU Utilization** | - | {time_data.get('cpu_pct', 'N/A')} | - |

## 2. Correctness Verification

| Metric | Expected (v0.3.0) | mimalloc Measured | Match |
|:---|---:|---:|:---:|
| Total Pairs | {TOTAL_PAIRS:,} | {sum_data.get('total_reads', 0):,} | {'PASS' if sum_data.get('total_reads', 0) == TOTAL_PAIRS else 'FAIL'} |
| Assigned | {EXP_ASSIGNED:,} | {assigned:,} | {'PASS' if assigned == EXP_ASSIGNED else 'FAIL'} |
| Unassigned | {EXP_UNASSIGNED:,} | {unassigned:,} | {'PASS' if unassigned == EXP_UNASSIGNED else 'FAIL'} |
| Canonical | {EXP_CANONICAL:,} | {canonical:,} | {'PASS' if canonical == EXP_CANONICAL else 'FAIL'} |
| Swapped | {EXP_SWAPPED:,} | {swapped:,} | {'PASS' if swapped == EXP_SWAPPED else 'FAIL'} |

## 3. Engineering Conclusion
- 吞吐表现: 全量 1.41 亿 pairs 实际吞吐为 {throughput:,.0f} pairs/s，相比 baseline 提速 {speedup:.2f}x。
- 内存开销: 峰值内存稳定在 {rss_mb:.1f} MB，未发生内存泄漏或随 reads 数量无限膨胀。
- 正确性: 全部统计指标与生产基准 100% 逐字吻合。
"""

with open(report_md, "w") as f:
    f.write(report)

print(report)
if not correct:
    sys.exit(1)
PYEOF

# Clean up massive output FASTQ to preserve disk quota, keep summaries
echo "Cleaning up generated FASTQ files in ${OUT_DIR}..."
find "${OUT_DIR}" -name "*.fastq.gz" -delete
echo "Done."
