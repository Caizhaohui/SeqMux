#!/usr/bin/env python3
"""v0.3 production thread tune, full I464/I395 runs, count compare, FASTQ checks.

Does not modify SeqMux. Production path:
  gzip input + default adapter trim + mismatch 0 + gzip compression level 1.
"""
import gzip
import os
import shutil
import statistics
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

BINARY = "target/release/seqmux"
TUNE_R1 = "tmp/upx_bench/i464_2000000_1.fq.gz"
TUNE_R2 = "tmp/upx_bench/i464_2000000_2.fq.gz"
TUNE_PAIRS = 2_000_000
TUNE_THREADS = [8, 12, 16]
REPS = 3

I464 = {
    "name": "I464",
    "r1": "/hpcfs/fhome/caizhh/03_PCR_analysis/25_464469erdaimix/E260702007_L01_464469erdaimix_1.fq.gz",
    "r2": "/hpcfs/fhome/caizhh/03_PCR_analysis/25_464469erdaimix/E260702007_L01_464469erdaimix_2.fq.gz",
    "barcodes": "tests/fixtures/I464-469erdai_barcode_and_name.csv",
    "python": "/hpcfs/fhome/caizhh/03_PCR_analysis/25_464469erdaimix/results/sample_read_counts.csv",
    "out": "tmp/v03_i464",
    "prefix": "i464",
    "expect_total": 140_771_720,
    "expect_assigned": 87_477_095,
    "expect_canonical": 45_393_500,
    "expect_swapped": 42_083_595,
    "bc1": "Barcode1",
}
I395 = {
    "name": "I395",
    "r1": "/hpcfs/fhome/caizhh/03_PCR_analysis/23_I395erdaimix/E260617005_L01_395erdaimix_1.fq.gz",
    "r2": "/hpcfs/fhome/caizhh/03_PCR_analysis/23_I395erdaimix/E260617005_L01_395erdaimix_2.fq.gz",
    "barcodes": "/hpcfs/fhome/caizhh/03_PCR_analysis/23_I395erdaimix/I395erdai_gst_barcode_and_name.csv",
    "python": "/hpcfs/fhome/caizhh/03_PCR_analysis/23_I395erdaimix/results/sample_read_counts.csv",
    "out": "tmp/v03_i395",
    "prefix": "i395",
    "expect_total": 60_523_000,
    "expect_assigned": 39_455_648,
    "expect_canonical": 18_606_371,
    "expect_swapped": 20_849_277,
}

# Must match src/cli.rs defaults used by the production command.
ADAPTER_R1 = b"AGATCGGAAGAGCACACGTCTGAA"
ADAPTER_R2 = b"AGATCGGAAGAGCGTCGTG"


def log(msg):
    print(msg, flush=True)


def check_env():
    log("=== ENVIRONMENT ===")
    log("hostname: " + subprocess.check_output(["hostname"], text=True).strip())
    log("nproc: " + subprocess.check_output(["nproc"], text=True).strip())
    log(subprocess.check_output(["taskset", "-cp", str(os.getpid())], text=True).strip())
    log("SLURM_CPUS_ON_NODE: " + os.environ.get("SLURM_CPUS_ON_NODE", "N/A"))
    log("====================")


def parse_time(stderr_text):
    wall = user = sys_s = rss = cpu = None
    for line in stderr_text.splitlines():
        line = line.strip()
        if "Elapsed (wall clock) time" in line:
            parts = line.split("):")[-1].strip().split(":")
            if len(parts) == 2:
                wall = float(parts[0]) * 60 + float(parts[1])
            elif len(parts) == 3:
                wall = float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])
            else:
                wall = float(parts[0])
        elif line.startswith("User time"):
            user = float(line.split(":")[-1].strip())
        elif line.startswith("System time"):
            sys_s = float(line.split(":")[-1].strip())
        elif "Percent of CPU" in line:
            cpu = float(line.split(":")[-1].strip().replace("%", ""))
        elif "Maximum resident set size" in line:
            rss = int(line.split(":")[-1].strip())
    return wall, user, sys_s, cpu, rss


def load_summary(path):
    data = {"samples": {}}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            k, v = parts[0], parts[1]
            if k.startswith("sample:"):
                data["samples"][k.split(":", 1)[1]] = int(v)
            elif k in {
                "total_reads",
                "assigned",
                "unassigned",
                "adapter_trimmed",
                "orientation_canonical",
                "orientation_swapped",
            }:
                data[k] = int(v)
    return data


def dir_size(path):
    total = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            fp = os.path.join(root, name)
            if not os.path.islink(fp):
                total += os.path.getsize(fp)
    return total


def run_demux(r1, r2, barcodes, out_dir, prefix, threads):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    cmd = [
        "/usr/bin/time",
        "-v",
        BINARY,
        "demux",
        "-i",
        r1,
        "-I",
        r2,
        "-b",
        barcodes,
        "-o",
        out_dir,
        "-p",
        prefix,
        "-t",
        str(threads),
        "--compression-level",
        "1",
        "--force",
    ]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        log("FAIL: " + " ".join(cmd))
        log(p.stderr[-3000:])
        sys.exit(1)
    wall, user, sys_s, cpu, rss = parse_time(p.stderr)
    summary = load_summary(os.path.join(out_dir, f"{prefix}.summary.tsv"))
    nbytes = dir_size(out_dir)
    return {
        "wall": wall,
        "user": user,
        "sys": sys_s,
        "cpu": cpu,
        "rss": rss,
        "assigned": summary.get("assigned"),
        "trimmed": summary.get("adapter_trimmed"),
        "total": summary.get("total_reads"),
        "unassigned": summary.get("unassigned"),
        "bytes": nbytes,
        "summary": summary,
    }


def median(vals):
    return statistics.median(vals)


def tune_threads():
    log("\n=== 1. PRODUCTION THREAD TUNING (2M I464, gzip in/out level 1) ===")
    rows = []
    for threads in TUNE_THREADS:
        runs = []
        for i in range(REPS):
            out = f"tmp/v03_tune/t{threads}_r{i}"
            m = run_demux(TUNE_R1, TUNE_R2, I464["barcodes"], out, "bench", threads)
            runs.append(m)
            log(
                f"  -t {threads} rep {i+1}: wall={m['wall']:.2f}s user={m['user']:.2f}s "
                f"sys={m['sys']:.2f}s cpu={m['cpu']:.0f}% thpt={TUNE_PAIRS/m['wall']:.0f} "
                f"rss={m['rss']} assigned={m['assigned']} trimmed={m['trimmed']} "
                f"bytes={m['bytes']}"
            )
            shutil.rmtree(out, ignore_errors=True)
        row = {
            "threads": threads,
            "wall": median([r["wall"] for r in runs]),
            "user": median([r["user"] for r in runs]),
            "sys": median([r["sys"] for r in runs]),
            "cpu": median([r["cpu"] for r in runs]),
            "thpt": median([TUNE_PAIRS / r["wall"] for r in runs]),
            "rss": median([r["rss"] for r in runs]),
            "assigned": runs[0]["assigned"],
            "trimmed": runs[0]["trimmed"],
            "bytes": median([r["bytes"] for r in runs]),
        }
        rows.append(row)
    log("\n| Threads | Wall | User | Sys | CPU % | pairs/s | Peak RSS | Assigned | Adapter trimmed | Output bytes |")
    log("| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |")
    for r in rows:
        log(
            f"| **-t {r['threads']}** | {r['wall']:.2f} s | {r['user']:.2f} s | {r['sys']:.2f} s | "
            f"{r['cpu']:.1f}% | {r['thpt']:,.0f} | {r['rss']:,.0f} kB | {r['assigned']:,} | "
            f"{r['trimmed']:,} | {r['bytes']:,.0f} |"
        )
    best = min(rows, key=lambda r: (r["wall"], r["threads"]))
    log(f"\nSELECTED_THREADS={best['threads']} (lowest median wall clock)")
    return best["threads"]


def compare_dataset(ds, summary_path):
    cmd = [
        "python3",
        "scripts/compare_python_counts.py",
        summary_path,
        ds["python"],
        "--expect-total",
        str(ds["expect_total"]),
        "--expect-assigned",
        str(ds["expect_assigned"]),
        "--expect-canonical",
        str(ds["expect_canonical"]),
        "--expect-swapped",
        str(ds["expect_swapped"]),
    ]
    log("\n--- count compare " + ds["name"] + " ---")
    p = subprocess.run(cmd, text=True)
    if p.returncode != 0:
        log(f"COUNT MISMATCH on {ds['name']}: investigate before release")
        sys.exit(2)


def load_barcodes(path):
    samples = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split(",")
        i_name = header.index("SampleNumber")
        i_b1 = header.index("Barcode1")
        i_b2 = header.index("Barcode2")
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split(",")
            samples[parts[i_name]] = (parts[i_b1].upper(), parts[i_b2].upper())
    return samples


def iter_fastq(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            if not qual:
                raise RuntimeError(f"truncated FASTQ {path}")
            yield header[1:].strip(), seq.strip(), qual.strip()


def pair_id(header):
    h = header.split()[0]
    if h.endswith("/1") or h.endswith("/2"):
        h = h[:-2]
    return h


def representative_samples(summary, k=3):
    items = sorted(summary["samples"].items(), key=lambda kv: kv[1])
    items = [(n, c) for n, c in items if c > 0]
    if not items:
        return []
    picks = [items[0], items[len(items) // 2], items[-1]]
    # unique preserve order
    seen = set()
    out = []
    for n, c in picks:
        if n not in seen:
            seen.add(n)
            out.append((n, c))
    return out[:k]


def validate_fastq(ds):
    log(f"\n=== 4. OUTPUT VALIDATION {ds['name']} ===")
    summary = load_summary(os.path.join(ds["out"], f"{ds['prefix']}.summary.tsv"))
    barcodes = load_barcodes(ds["barcodes"])
    samples = representative_samples(summary)
    log("representative samples: " + ", ".join(f"{n}={c}" for n, c in samples))

    watch = {}  # pair_id -> sample
    file_info = {}
    for name, count in samples:
        r1 = os.path.join(ds["out"], f"{ds['prefix']}_{name}_R1.fastq.gz")
        r2 = os.path.join(ds["out"], f"{ds['prefix']}_{name}_R2.fastq.gz")
        for path in (r1, r2):
            p = subprocess.run(["gzip", "-t", path], capture_output=True, text=True)
            if p.returncode != 0:
                log(f"gzip -t FAIL {path}: {p.stderr}")
                sys.exit(3)
        n1 = n2 = 0
        watched_out = []
        first_pairs = []
        with gzip.open(r1, "rt") as f1, gzip.open(r2, "rt") as f2:
            while True:
                h1 = f1.readline()
                if not h1:
                    break
                s1 = f1.readline().strip()
                f1.readline()
                f1.readline()
                h2 = f2.readline()
                s2 = f2.readline().strip()
                f2.readline()
                f2.readline()
                if not h2:
                    log(f"R2 ended early for {name}")
                    sys.exit(3)
                id1, id2 = pair_id(h1[1:].strip()), pair_id(h2[1:].strip())
                if id1 != id2:
                    log(f"ID mismatch {name}: {id1} vs {id2}")
                    sys.exit(3)
                n1 += 1
                n2 += 1
                # Head and tail of the output file must be a subsequence of input order.
                if n1 <= 30 or n1 > count - 20:
                    if id1 not in watch:
                        watch[id1] = name
                        watched_out.append(id1)
                    if n1 <= 30:
                        first_pairs.append((id1, s1, s2))
        if n1 != count or n2 != count:
            log(f"record count {name}: files={n1}/{n2} summary={count}")
            sys.exit(3)
        # checksums
        c1 = subprocess.check_output(["sha256sum", r1], text=True).split()[0]
        c2 = subprocess.check_output(["sha256sum", r2], text=True).split()[0]
        file_info[name] = {
            "watched_out": watched_out,
            "first": first_pairs,
            "sha_r1": c1,
            "sha_r2": c2,
            "n": n1,
        }
        log(f"  {name}: paired records={n1} gzip_ok sha256_R1={c1[:16]}… sha256_R2={c2[:16]}…")

    # One pass over inputs: confirm watched IDs keep output order and canonical trim.
    order_pos = {name: [] for name, _ in samples}
    found = {}
    log("  scanning inputs for order + orientation spot checks...")
    r1_iter = iter_fastq(ds["r1"])
    r2_iter = iter_fastq(ds["r2"])
    pos = 0
    remaining = set(watch)
    while remaining:
        try:
            h1, s1, _q1 = next(r1_iter)
            h2, s2, _q2 = next(r2_iter)
        except StopIteration:
            break
        pid = pair_id(h1)
        if pair_id(h2) != pid:
            log(f"input mates disagree at {pid}")
            sys.exit(3)
        pos += 1
        if pid in remaining:
            name = watch[pid]
            order_pos[name].append(pid)
            found[pid] = (s1, s2)
            remaining.discard(pid)
        if pos % 20_000_000 == 0:
            log(f"    scanned {pos:,} pairs, watching {len(remaining)}")
    if remaining:
        log(f"missing watched IDs: {list(remaining)[:5]}")
        sys.exit(3)

    adapter_suffix_ok = 0
    adapter_suffix_n = 0
    for name, _count in samples:
        watched_out = file_info[name]["watched_out"]
        if order_pos[name] != watched_out:
            log(f"ORDER FAIL {name}: output ID order != input appearance order")
            log(f"  output head: {watched_out[:3]} tail: {watched_out[-3:]}")
            log(f"  input  head: {order_pos[name][:3]} tail: {order_pos[name][-3:]}")
            sys.exit(3)
        bc1, bc2 = barcodes[name]
        for pid, out_r1, out_r2 in file_info[name]["first"]:
            in_r1, in_r2 = found[pid]
            if in_r1.startswith(bc1) and in_r2.startswith(bc2):
                src_r1, src_r2 = in_r1[len(bc1) :], in_r2[len(bc2) :]
            elif in_r1.startswith(bc2) and in_r2.startswith(bc1):
                src_r1, src_r2 = in_r2[len(bc1) :], in_r1[len(bc2) :]
            else:
                log(f"ORIENTATION FAIL {name} {pid}: barcodes not at 5'")
                sys.exit(3)
            if not src_r1.startswith(out_r1) or not src_r2.startswith(out_r2):
                log(f"TRIM FAIL {name} {pid}: output is not a prefix of barcode-stripped canonical mates")
                sys.exit(3)
            for src, out, adapter in (
                (src_r1, out_r1, ADAPTER_R1),
                (src_r2, out_r2, ADAPTER_R2),
            ):
                if len(out) < len(src):
                    removed = src[len(out) :].upper()
                    adapter_suffix_n += 1
                    if adapter.startswith(removed) or removed.startswith(adapter[: min(8, len(removed))]):
                        adapter_suffix_ok += 1
        log(
            f"  {name}: order preserved ({len(watched_out)} head/tail IDs), canonical prefix OK"
        )

    # Tail order: collect last 20 IDs and confirm they appear in that order in a second
    # observation already stored only if they were watched. For large files last IDs
    # were not watched. Scan is done. Verify last-20 IDs are a subsequence by a focused
    # second pass only if the sample is small enough that we included them in watch.
    log(
        f"  3' suffix spot-check: {adapter_suffix_ok}/{adapter_suffix_n} removed tails match adapter prefix"
    )
    log(f"PASS output validation {ds['name']}")


def full_run(ds, threads):
    log(f"\n=== FULL {ds['name']} production run -t {threads} ===")
    m = run_demux(ds["r1"], ds["r2"], ds["barcodes"], ds["out"], ds["prefix"], threads)
    pairs = m["total"] or 0
    log(
        f"total_pairs={pairs:,} wall={m['wall']:.2f}s user={m['user']:.2f}s sys={m['sys']:.2f}s "
        f"cpu={m['cpu']:.1f}% thpt={pairs / m['wall']:.0f} rss={m['rss']:,} kB "
        f"assigned={m['assigned']:,} unassigned={m['unassigned']:,} "
        f"adapter_trimmed={m['trimmed']:,} output_bytes={m['bytes']:,}"
    )
    log("per-sample counts:")
    for name, count in sorted(m["summary"]["samples"].items()):
        log(f"  {name}\t{count}")
    summary_path = os.path.join(ds["out"], f"{ds['prefix']}.summary.tsv")
    sha = subprocess.check_output(["sha256sum", summary_path], text=True).strip()
    log("summary_sha256: " + sha)
    compare_dataset(ds, summary_path)
    validate_fastq(ds)
    return m


def main():
    check_env()
    if not os.path.exists(BINARY):
        log("missing " + BINARY)
        sys.exit(1)
    threads = tune_threads()
    with open("tmp/v03_selected_threads.txt", "w") as f:
        f.write(str(threads) + "\n")
    full_run(I464, threads)
    full_run(I395, threads)
    log("\nV03_VALIDATION_OK")


if __name__ == "__main__":
    main()
