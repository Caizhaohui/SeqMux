#!/bin/bash
#SBATCH --job-name=seqmux-m17a
#SBATCH --partition=qcpu_18i
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=benchmarks/m17a_ab_benchmark/slurm_%j.out
#SBATCH --error=benchmarks/m17a_ab_benchmark/slurm_%j.err

set -euo pipefail

if [[ -z ${SLURM_JOB_ID:-} ]]; then
    echo "This script must be submitted via sbatch." >&2
    exit 2
fi

cd "${SLURM_SUBMIT_DIR:?Submit from the SeqMux repository root}"

echo "=================================================="
echo "SeqMux M17A: Strict A/B Release Profile & Allocator Benchmark"
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $SLURM_NODELIST"
echo "CPUs:   $SLURM_CPUS_PER_TASK"
echo "Date:   $(date -Iseconds)"
echo "=================================================="

mkdir -p benchmarks/m17a_ab_benchmark

python3 scripts/run_m17a_ab_benchmark.py

echo ""
echo "=================================================="
echo "M17A Benchmark Completed"
echo "Date: $(date -Iseconds)"
echo "=================================================="
