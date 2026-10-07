#!/bin/bash
#SBATCH --job-name=seqmux-m17
#SBATCH --partition=qcpu_18i
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-4%2

set -euo pipefail

if [[ -z ${SLURM_JOB_ID:-} || -z ${SLURM_ARRAY_TASK_ID:-} ]]; then
    echo "Run through sbatch with array tasks 0-4; full datasets require a compute allocation." >&2
    exit 2
fi

binary=${1:?Provide an absolute path to the validated M17 release binary}
run_root=${2:?Provide an absolute path to a fresh output directory}
[[ $binary = /* && $run_root = /* && -x $binary ]]
cd "${SLURM_SUBMIT_DIR:?Submit from the SeqMux repository root}"

dataset=i464
mismatches=0
delta=0
case "$SLURM_ARRAY_TASK_ID" in
    0) label=i464_exact ;;
    1) label=i395_exact; dataset=i395 ;;
    2|3|4) delta=$((SLURM_ARRAY_TASK_ID - 2)); mismatches=1; label=i464_mm1_delta${delta} ;;
    *) echo "Expected array task 0-4" >&2; exit 2 ;;
esac

if [[ $dataset = i464 ]]; then
    data=/hpcfs/fhome/caizhh/03_PCR_analysis/25_464469erdaimix
    r1=$data/E260702007_L01_464469erdaimix_1.fq.gz
    r2=$data/E260702007_L01_464469erdaimix_2.fq.gz
    barcodes=tests/fixtures/I464-469erdai_barcode_and_name.csv
    total=140771720; assigned=87477095; unassigned=53294625
    canonical=45393500; swapped=42083595
else
    data=/hpcfs/fhome/caizhh/03_PCR_analysis/23_I395erdaimix
    r1=$data/E260617005_L01_395erdaimix_1.fq.gz
    r2=$data/E260617005_L01_395erdaimix_2.fq.gz
    barcodes=$data/I395erdai_gst_barcode_and_name.csv
    total=60523000; assigned=39455648; unassigned=21067352
    canonical=18606371; swapped=20849277
fi

out=$run_root/$label
mkdir -p "$run_root"
mkdir "$out"
args=(demux -i "$r1" -I "$r2" -b "$barcodes" -o "$out" -p "$dataset"
    -t 12 --orientation both --compression-level 1 -q 0
    --mismatches-1 "$mismatches" --mismatches-2 "$mismatches")
if [[ $mismatches = 1 ]]; then
    args+=(--counts-only --min-mismatch-delta "$delta")
fi

hostname > "$out/node.txt"
sha256sum "$binary" > "$out/binary.sha256"
printf '%q ' "$binary" "${args[@]}" > "$out/command.txt"
printf '\n' >> "$out/command.txt"
/usr/bin/time -v "$binary" "${args[@]}" > "$out/stdout.log" 2> "$out/stderr.log"
summary=$out/$dataset.summary.tsv

if [[ $mismatches = 0 ]]; then
    python3 scripts/compare_python_counts.py "$summary" "$data/results/sample_read_counts.csv" \
        --expect-total "$total" --expect-assigned "$assigned" \
        --expect-canonical "$canonical" --expect-swapped "$swapped" > "$out/python_parity.txt"
    awk -F '\t' -v t="$total" -v a="$assigned" -v u="$unassigned" -v c="$canonical" -v s="$swapped" '
        { counts[$1] = $2 }
        END {
            if (counts["total_reads"] != t || counts["assigned"] != a ||
                counts["unassigned"] != u || counts["orientation_canonical"] != c ||
                counts["orientation_swapped"] != s || counts["ambiguous"] != 0) exit 1
        }
    ' "$summary"
    baseline=tmp/v03_$dataset/$dataset.summary.tsv
elif [[ $delta = 0 ]]; then
    baseline=tmp/m13_i464_mm1/i464mm1.summary.tsv
else
    baseline=
fi

if [[ -n $baseline ]]; then
    awk -F '\t' '
        function tracked(k) {
            return k ~ /^sample:/ || k ~ /^(total_reads|assigned|unassigned|ambiguous|orientation_canonical|orientation_swapped)$/
        }
        NR == FNR { if (tracked($1)) old[$1] = $2; next }
        tracked($1) {
            seen[$1] = 1
            if (!($1 in old) || old[$1] != $2) { print "BASELINE MISMATCH", $1, old[$1], $2; bad = 1 }
        }
        END {
            for (k in old) if (!(k in seen)) { print "MISSING METRIC", k; bad = 1 }
            exit bad
        }
    ' "$baseline" "$summary" > "$out/baseline_comparison.txt"
fi

printf 'PASS %s\n' "$label" | tee "$out/PASS"
