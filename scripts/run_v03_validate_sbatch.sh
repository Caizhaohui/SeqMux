#!/bin/bash
#SBATCH --job-name=seqmux-v03
#SBATCH --partition=qcpu_18i
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/03_SeqMux/tmp/slurm_v03-%j.out
#SBATCH --error=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/03_SeqMux/tmp/slurm_v03-%j.err

set -euo pipefail
cd /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/03_SeqMux

echo "Running on node: $(hostname)"
echo "Available CPUs: $(nproc)"
echo "CPU affinity: $(taskset -cp $$)"

cargo build --release -q
PYTHONUNBUFFERED=1 python3 -u scripts/v03_production_validate.py
