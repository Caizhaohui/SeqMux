#!/usr/bin/env python3
"""SeqMux M17A — Release Build / Allocator Strict A/B Benchmark.

Evaluates 4 runtime variants:
  A — Current baseline (thin LTO, codegen-units=16, default allocator)
  B — Fat LTO only (fat LTO, codegen-units=1, default allocator)
  C — mimalloc only (thin LTO, codegen-units=16, mimalloc)
  D — Fat LTO + mimalloc (fat LTO, codegen-units=1, mimalloc)

Follows strict interleaved benchmark methodology across 5 measured rounds (+1 warmup)
on the primary I464 2,000,000 paired-read production workload (-t 12, --compression-level 1),
plus diagnostic counts-only workload.
"""

import hashlib
import json
import os
import shutil
import socket
import statistics
import subprocess
import sys
import time

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BENCH_DIR = os.path.join(ROOT_DIR, "benchmarks/m17a_ab_benchmark")
PROV_DIR = os.path.join(BENCH_DIR, "provenance")
BIN_DIR = os.path.join(BENCH_DIR, "binaries")
RUNS_DIR = os.path.join(ROOT_DIR, "tmp/m17a_runs")

I464_R1 = os.path.join(ROOT_DIR, "tmp/upx_bench/i464_2000000_1.fq.gz")
I464_R2 = os.path.join(ROOT_DIR, "tmp/upx_bench/i464_2000000_2.fq.gz")
I464_BARCODES = os.path.join(ROOT_DIR, "tests/fixtures/I464-469erdai_barcode_and_name.csv")

EXPECTED_COUNTS = {
    "total_reads": 2_000_000,
    "assigned": 1_251_876,
    "unassigned": 748_124,
    "orientation_canonical": 629_825,
    "orientation_swapped": 622_051,
    "adapter_trimmed": 16_542,
}

VARIANTS = {
    "A": {
        "name": "A_baseline",
        "description": "Current baseline (thin LTO, codegen-units=16, default allocator)",
        "cargo_args": ["build", "--release", "--locked"],
        "target_bin": os.path.join(ROOT_DIR, "target/release/seqmux"),
    },
    "B": {
        "name": "B_fat_lto",
        "description": "Fat LTO only (fat LTO, codegen-units=1, default allocator)",
        "cargo_args": ["build", "--profile", "release-fat", "--locked"],
        "target_bin": os.path.join(ROOT_DIR, "target/release-fat/seqmux"),
    },
    "C": {
        "name": "C_mimalloc",
        "description": "mimalloc only (thin LTO, codegen-units=16, mimalloc)",
        "cargo_args": ["build", "--release", "--features", "mimalloc-allocator", "--locked"],
        "target_bin": os.path.join(ROOT_DIR, "target/release/seqmux"),
    },
    "D": {
        "name": "D_both",
        "description": "Fat LTO + mimalloc (fat LTO, codegen-units=1, mimalloc)",
        "cargo_args": ["build", "--profile", "release-fat", "--features", "mimalloc-allocator", "--locked"],
        "target_bin": os.path.join(ROOT_DIR, "target/release-fat/seqmux"),
    },
}

ROUNDS_ORDER = [
    ["A", "B", "C", "D"],  # Round 1
    ["D", "C", "B", "A"],  # Round 2
    ["B", "D", "A", "C"],  # Round 3
    ["C", "A", "D", "B"],  # Round 4
    ["A", "C", "B", "D"],  # Round 5
]


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def parse_time_output(time_str):
    metrics = {
        "wall_seconds": 0.0,
        "user_seconds": 0.0,
        "sys_seconds": 0.0,
        "cpu_percent": 0.0,
        "peak_rss_kb": 0,
        "exit_code": 0,
    }
    for line in time_str.splitlines():
        line = line.strip()
        if "User time (seconds):" in line:
            metrics["user_seconds"] = float(line.split(":")[-1].strip())
        elif "System time (seconds):" in line:
            metrics["sys_seconds"] = float(line.split(":")[-1].strip())
        elif "Percent of CPU this job got:" in line:
            metrics["cpu_percent"] = float(line.split(":")[-1].strip().rstrip("%"))
        elif "Maximum resident set size (kbytes):" in line:
            metrics["peak_rss_kb"] = int(line.split(":")[-1].strip())
        elif "Exit status:" in line:
            metrics["exit_code"] = int(line.split(":")[-1].strip())
        elif "Elapsed (wall clock) time (h:mm:ss or m:ss):" in line:
            raw = line.split("):")[-1].strip()
            parts = raw.split(":")
            if len(parts) == 3:
                metrics["wall_seconds"] = float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])
            elif len(parts) == 2:
                metrics["wall_seconds"] = float(parts[0]) * 60 + float(parts[1])
            else:
                metrics["wall_seconds"] = float(parts[0])
    return metrics


def parse_summary_tsv(summary_path):
    metrics = {}
    if not os.path.exists(summary_path):
        return metrics
    with open(summary_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                k, v = parts[0], parts[1]
                try:
                    metrics[k] = int(v) if "." not in v else float(v)
                except ValueError:
                    metrics[k] = v
    return metrics


def freeze_reference():
    os.makedirs(PROV_DIR, exist_ok=True)
    log("Freezing reference provenance...")

    def run_cmd_to_file(cmd, out_file):
        try:
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            with open(os.path.join(PROV_DIR, out_file), "w") as f:
                f.write(res.stdout)
        except Exception as e:
            with open(os.path.join(PROV_DIR, out_file), "w") as f:
                f.write(f"Error: {e}\n")

    run_cmd_to_file(["git", "rev-parse", "HEAD"], "git_commit.txt")
    run_cmd_to_file(["git", "status", "--porcelain"], "git_status.txt")
    run_cmd_to_file(["rustc", "-Vv"], "rustc_version.txt")
    run_cmd_to_file(["cargo", "-V"], "cargo_version.txt")
    with open(os.path.join(PROV_DIR, "cargo_lock_sha256.txt"), "w") as f:
        f.write(sha256_file(os.path.join(ROOT_DIR, "Cargo.lock")) + "  Cargo.lock\n")
    with open(os.path.join(PROV_DIR, "hostname.txt"), "w") as f:
        f.write(socket.gethostname() + "\n")
    run_cmd_to_file(["lscpu"], "lscpu.txt")
    run_cmd_to_file(["nproc"], "nproc.txt")

    # CPU affinity / Cpus_allowed
    try:
        with open("/proc/self/status") as f:
            lines = [l for l in f if "Cpus_allowed" in l]
        with open(os.path.join(PROV_DIR, "cpus_allowed.txt"), "w") as f:
            f.writelines(lines)
    except Exception:
        pass

    # SLURM env
    with open(os.path.join(PROV_DIR, "slurm_env.txt"), "w") as f:
        for k, v in os.environ.items():
            if k.startswith("SLURM_"):
                f.write(f"{k}={v}\n")


def build_variants():
    os.makedirs(BIN_DIR, exist_ok=True)
    build_metrics = {}

    log("Building 4 runtime variants...")
    for var_id, var in VARIANTS.items():
        vname = var["name"]
        dst_bin = os.path.join(BIN_DIR, f"seqmux_{vname}")
        log(f"--- Building Variant {var_id} ({vname}): {var['description']} ---")

        # Clean seqmux crate to measure clean compile time
        subprocess.run(["cargo", "clean", "-p", "seqmux"], cwd=ROOT_DIR, check=True)

        t0 = time.perf_counter()
        res = subprocess.run(["cargo"] + var["cargo_args"], cwd=ROOT_DIR, capture_output=True, text=True)
        t1 = time.perf_counter()

        if res.returncode != 0:
            log(f"ERROR: Build failed for {vname}:\n{res.stderr}")
            sys.exit(1)

        build_time = t1 - t0
        src_bin = var["target_bin"]
        shutil.copy2(src_bin, dst_bin)

        bin_size = os.path.getsize(dst_bin)
        bin_sha = sha256_file(dst_bin)

        build_metrics[var_id] = {
            "name": vname,
            "description": var["description"],
            "build_time_seconds": build_time,
            "binary_size_bytes": bin_size,
            "sha256": bin_sha,
            "binary_path": dst_bin,
        }
        log(f"Variant {var_id} built in {build_time:.2f}s, size: {bin_size / 1024 / 1024:.2f} MB, sha256: {bin_sha[:16]}...")

    with open(os.path.join(BENCH_DIR, "build_metrics.json"), "w") as f:
        json.dump(build_metrics, f, indent=2)
    return build_metrics


def run_benchmark_workload(build_metrics, is_diagnostic=False):
    workload_name = "diagnostic_counts_only" if is_diagnostic else "primary_production"
    log(f"==================================================")
    log(f"Starting {workload_name.upper()} benchmark")
    log(f"==================================================")

    # 1. Warm-up round
    log("--- Running 1 Warm-up per variant ---")
    for var_id in ["A", "B", "C", "D"]:
        bin_path = build_metrics[var_id]["binary_path"]
        out_dir = os.path.join(RUNS_DIR, f"warmup_{workload_name}_{var_id}")
        shutil.rmtree(out_dir, ignore_errors=True)
        cmd = [
            bin_path, "demux",
            "-i", I464_R1, "-I", I464_R2,
            "-b", I464_BARCODES,
            "-o", out_dir, "-p", "test",
            "-t", "12",
            "--force",
        ]
        if is_diagnostic:
            cmd.append("--counts-only")
        else:
            cmd.extend(["--compression-level", "1"])

        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        shutil.rmtree(out_dir, ignore_errors=True)
        log(f"Warm-up {var_id} completed.")

    # 2. Five Interleaved Measured Rounds
    results = {v: [] for v in ["A", "B", "C", "D"]}

    for round_idx, round_vars in enumerate(ROUNDS_ORDER, 1):
        log(f"\n--- Round {round_idx}/5: Order {' -> '.join(round_vars)} ---")
        for var_id in round_vars:
            bin_path = build_metrics[var_id]["binary_path"]
            out_dir = os.path.join(RUNS_DIR, f"{workload_name}_{var_id}_rep{round_idx}")
            shutil.rmtree(out_dir, ignore_errors=True)
            time_log = out_dir + ".time"

            cmd = [
                "/usr/bin/time", "-v", "-o", time_log,
                bin_path, "demux",
                "-i", I464_R1, "-I", I464_R2,
                "-b", I464_BARCODES,
                "-o", out_dir, "-p", "test",
                "-t", "12",
                "--force",
            ]
            if is_diagnostic:
                cmd.append("--counts-only")
            else:
                cmd.extend(["--compression-level", "1"])

            t0 = time.perf_counter()
            p = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            t1 = time.perf_counter()

            if p.returncode != 0:
                log(f"ERROR: Run failed for Variant {var_id} in round {round_idx}")
                sys.exit(1)

            with open(time_log) as f:
                metrics = parse_time_output(f.read())
            if metrics["wall_seconds"] == 0.0:
                metrics["wall_seconds"] = t1 - t0

            # Output size
            total_bytes = 0
            if os.path.exists(out_dir):
                for root, _, files in os.walk(out_dir):
                    for fname in files:
                        total_bytes += os.path.getsize(os.path.join(root, fname))
            metrics["output_bytes"] = total_bytes

            # Parse summary
            summary_file = os.path.join(out_dir, "test.summary.tsv")
            summary = parse_summary_tsv(summary_file)
            metrics["summary"] = summary

            # Correctness assertion
            for k, expected_v in EXPECTED_COUNTS.items():
                actual_v = summary.get(k)
                if actual_v != expected_v:
                    log(f"FATAL: Variant {var_id} count mismatch on {k}: got {actual_v}, expected {expected_v}")
                    sys.exit(2)

            pairs_s = EXPECTED_COUNTS["total_reads"] / metrics["wall_seconds"]
            metrics["pairs_per_second"] = pairs_s

            results[var_id].append(metrics)
            log(f"  Variant {var_id}: wall={metrics['wall_seconds']:.2f}s | {pairs_s/1000:.1f}k pairs/s | rss={metrics['peak_rss_kb']/1024:.1f} MB (PASS)")

            # For production workload, preserve rep1 for FASTQ decompression parity check, delete others
            if not (not is_diagnostic and round_idx == 1):
                shutil.rmtree(out_dir, ignore_errors=True)
            if os.path.exists(time_log):
                os.remove(time_log)

    return results


def check_decompressed_parity():
    log("\n==================================================")
    log("Correctness Gate: Decompressed FASTQ Parity Check")
    log("==================================================")
    base_dir = os.path.join(RUNS_DIR, "primary_production_A_rep1")

    for var_id in ["B", "C", "D"]:
        cand_dir = os.path.join(RUNS_DIR, f"primary_production_{var_id}_rep1")
        log(f"Comparing Variant {var_id} against Baseline A...")

        base_files = sorted([f for f in os.listdir(base_dir) if f.endswith(".fastq.gz")])
        cand_files = sorted([f for f in os.listdir(cand_dir) if f.endswith(".fastq.gz")])

        if base_files != cand_files:
            log(f"FATAL: Output file lists differ between A and {var_id}!")
            sys.exit(3)

        for fname in base_files:
            base_gz = os.path.join(base_dir, fname)
            cand_gz = os.path.join(cand_dir, fname)

            # Compute decompressed sha256
            cmd_base = f"gzip -dc '{base_gz}' | sha256sum"
            cmd_cand = f"gzip -dc '{cand_gz}' | sha256sum"

            sha_base = subprocess.check_output(cmd_base, shell=True, text=True).split()[0]
            sha_cand = subprocess.check_output(cmd_cand, shell=True, text=True).split()[0]

            if sha_base != sha_cand:
                log(f"FATAL: Decompressed FASTQ mismatch on {fname} for Variant {var_id}!")
                sys.exit(4)

        log(f"Variant {var_id}: All {len(base_files)} sample FASTQ streams matched Baseline A byte-for-byte!")

    # Cleanup preserved rep1 directories
    for var_id in ["A", "B", "C", "D"]:
        shutil.rmtree(os.path.join(RUNS_DIR, f"primary_production_{var_id}_rep1"), ignore_errors=True)
    log("Decompressed FASTQ parity verification: ALL CHECKS PASSED.")


def generate_report(build_metrics, prod_results, diag_results):
    log("\nGenerating M17A Benchmark Report...")

    # Calculate statistics
    stats_table = []
    base_prod_wall = statistics.median([m["wall_seconds"] for m in prod_results["A"]])
    base_diag_wall = statistics.median([m["wall_seconds"] for m in diag_results["A"]])

    for var_id in ["A", "B", "C", "D"]:
        p_walls = [m["wall_seconds"] for m in prod_results[var_id]]
        p_pairs = [m["pairs_per_second"] for m in prod_results[var_id]]
        p_user = [m["user_seconds"] for m in prod_results[var_id]]
        p_sys = [m["sys_seconds"] for m in prod_results[var_id]]
        p_rss = [m["peak_rss_kb"] / 1024 for m in prod_results[var_id]]

        d_walls = [m["wall_seconds"] for m in diag_results[var_id]]

        med_p_wall = statistics.median(p_walls)
        med_p_pairs = statistics.median(p_pairs)
        speedup_prod = base_prod_wall / med_p_wall

        med_d_wall = statistics.median(d_walls)
        speedup_diag = base_diag_wall / med_d_wall

        stats_table.append({
            "variant": var_id,
            "name": VARIANTS[var_id]["name"],
            "prod_median_wall": med_p_wall,
            "prod_min_max": f"{min(p_walls):.2f}–{max(p_walls):.2f}",
            "prod_pairs_s": med_p_pairs,
            "prod_speedup": speedup_prod,
            "prod_user": statistics.median(p_user),
            "prod_sys": statistics.median(p_sys),
            "prod_rss": statistics.median(p_rss),
            "diag_median_wall": med_d_wall,
            "diag_speedup": speedup_diag,
            "bin_size_mb": build_metrics[var_id]["binary_size_bytes"] / 1024 / 1024,
            "build_time_s": build_metrics[var_id]["build_time_seconds"],
        })

    # Decision rules
    speedup_b = stats_table[1]["prod_speedup"]
    speedup_c = stats_table[2]["prod_speedup"]
    speedup_d = stats_table[3]["prod_speedup"]

    rss_a = stats_table[0]["prod_rss"]
    rss_c = stats_table[2]["prod_rss"]

    fat_lto_verdict = "KEEP" if speedup_b >= 1.02 else "REJECT"
    mimalloc_verdict = "KEEP" if (speedup_c >= 1.02 and rss_c <= rss_a * 1.5) else "REJECT"

    best_cfg = "A (Baseline thin LTO)"
    if fat_lto_verdict == "KEEP" and mimalloc_verdict == "KEEP":
        best_cfg = "D (Fat LTO + mimalloc)" if speedup_d >= max(speedup_b, speedup_c) else ("B" if speedup_b >= speedup_c else "C")
    elif fat_lto_verdict == "KEEP":
        best_cfg = "B (Fat LTO)"
    elif mimalloc_verdict == "KEEP":
        best_cfg = "C (mimalloc)"

    report_lines = [
        "# SeqMux M17A Benchmark Report — Release Build / Allocator Strict A/B",
        "",
        "## 1. Environment & Hardware",
        "",
        f"- **Date**: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"- **Host**: `{socket.gethostname()}`",
        f"- **SLURM Job**: `{os.environ.get('SLURM_JOB_ID', 'N/A')}`",
        f"- **CPUs allocated**: `{os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count())}`",
        f"- **Dataset**: I464 2,000,000 paired-read benchmark dataset (`tmp/upx_bench/i464_2000000_*.fq.gz`)",
        f"- **Command**: `-t 12 --compression-level 1` (production) / `-t 12 --counts-only` (diagnostic)",
        f"- **Methodology**: 1 warm-up + 5 interleaved measured repetitions across 4 variants",
        "",
        "## 2. Production Workload Results (Median of 5 Replicates)",
        "",
        "| Variant | Description | Wall (s) | Min–Max (s) | pairs/s | Speedup | User CPU | Sys CPU | Peak RSS (MB) |",
        "|:---|:---|---:|---:|---:|---:|---:|---:|---:|",
    ]

    for row in stats_table:
        report_lines.append(
            f"| **{row['variant']}** | {row['name']} | {row['prod_median_wall']:.2f} | {row['prod_min_max']} | {row['prod_pairs_s']:,.0f} | **{row['prod_speedup']:.2f}x** | {row['prod_user']:.2f}s | {row['prod_sys']:.2f}s | {row['prod_rss']:.1f} |"
        )

    report_lines.extend([
        "",
        "## 3. Diagnostic Counts-Only Workload & Build Trade-Offs",
        "",
        "| Variant | Diagnostic Wall (s) | Diag Speedup | Clean Build Time (s) | Binary Size (MB) |",
        "|:---|---:|---:|---:|---:|",
    ])

    for row in stats_table:
        report_lines.append(
            f"| **{row['variant']}** | {row['diag_median_wall']:.2f} | {row['diag_speedup']:.2f}x | {row['build_time_s']:.1f}s | {row['bin_size_mb']:.2f} MB |"
        )

    report_lines.extend([
        "",
        "## 4. Correctness Gate Verification",
        "",
        "- **Sample Counts**: All 4 variants reproduced exactly 1,251,876 assigned, 748,124 unassigned, 629,825 canonical, 622,051 swapped, and 16,542 adapter-trimmed pairs across all 35 samples.",
        "- **Decompressed FASTQ Content**: Verified byte-identical SHA-256 for all 70 sample FASTQ streams between Variant A and Variants B, C, D.",
        "- **Status**: **ALL CORRECTNESS CHECKS PASSED**.",
        "",
        "## 5. M17A Final Verdict",
        "",
        "```text",
        f"fat LTO:              {fat_lto_verdict} (speedup: {speedup_b:.2f}x, build time: {stats_table[1]['build_time_s']:.1f}s vs {stats_table[0]['build_time_s']:.1f}s)",
        f"mimalloc:             {mimalloc_verdict} (speedup: {speedup_c:.2f}x, RSS: {stats_table[2]['prod_rss']:.1f} MB vs {stats_table[0]['prod_rss']:.1f} MB)",
        f"RECOMMENDED BUILD:    {best_cfg}",
        "```",
        "",
    ])

    report_content = "\n".join(report_lines)
    report_path = os.path.join(BENCH_DIR, "M17A_REPORT.md")
    with open(report_path, "w") as f:
        f.write(report_content)

    print("\n" + "=" * 80)
    print(report_content)
    print("=" * 80)
    log(f"Report saved to {report_path}")


def main():
    freeze_reference()
    build_metrics = build_variants()
    prod_results = run_benchmark_workload(build_metrics, is_diagnostic=False)
    diag_results = run_benchmark_workload(build_metrics, is_diagnostic=True)
    check_decompressed_parity()
    generate_report(build_metrics, prod_results, diag_results)


if __name__ == "__main__":
    main()
