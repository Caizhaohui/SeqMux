#!/usr/bin/env python3
"""Compare SeqMux summary.tsv sample counts to Ultraplex output FASTQ files."""
from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path


def load_seqmux(path: Path):
    total = assigned = unassigned = None
    samples: dict[str, int] = {}
    with path.open() as f:
        next(f)
        for line in f:
            k, v = line.rstrip("\n").split("\t")
            if k in {
                "match_rate",
                "assigned_pct",
            } or k.endswith("_rate") or k.endswith("_pct"):
                continue
            try:
                n = int(v)
            except ValueError:
                continue
            if k == "total_reads":
                total = n
            elif k == "assigned":
                assigned = n
            elif k == "unassigned":
                unassigned = n
            elif k.startswith("sample:"):
                samples[k.split(":", 1)[1]] = n
    return total, assigned, unassigned, samples


def count_fastq_records(path: Path) -> int:
    opener = gzip.open if path.name.endswith(".gz") else open
    n = 0
    with opener(path, "rt") as f:
        for i, _ in enumerate(f, 1):
            pass
        n = i if i else 0
    if n % 4 != 0:
        raise SystemExit(f"{path}: {n} lines not divisible by 4")
    return n // 4


def load_ultraplex(out_dir: Path, prefix: str) -> dict[str, int]:
    """Count pairs from renamed `{prefix}_{sample}_Fwd.fastq.gz` files.

    Also accepts unrenamed `ultraplex_{prefix}_5bc_*_3bc_*_Fwd.fastq.gz` when
    sample names were not applied; those keys stay as barcode tokens.
    """
    samples: dict[str, int] = {}
    for path in sorted(out_dir.iterdir()):
        name = path.name
        if "tmp_thread" in name or name.endswith(".log"):
            continue
        for ext in (".fastq.gz", ".fastq"):
            if not name.endswith(ext):
                continue
            stem = name[: -len(ext)]
            # Renamed: upx_Sample_Fwd
            if stem.startswith(f"{prefix}_") and stem.endswith("_Fwd"):
                sample = stem[len(prefix) + 1 : -len("_Fwd")]
                if "no_match" in sample:
                    continue
                samples[sample] = count_fastq_records(path)
                break
            # Unrenamed: ultraplex_upx_5bc_XXX_3bc_YYY_Fwd
            marker = f"ultraplex_{prefix}_"
            if stem.startswith(marker) and stem.endswith("_Fwd"):
                sample = stem[len(marker) : -len("_Fwd")]
                if "no_match" in sample:
                    continue
                samples[sample] = count_fastq_records(path)
                break
    return samples


def collapse_swapped(samples: dict[str, int]) -> dict[str, int]:
    collapsed: dict[str, int] = defaultdict(int)
    for name, n in samples.items():
        base = name[: -len("_swapped")] if name.endswith("_swapped") else name
        collapsed[base] += n
    return dict(collapsed)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("seqmux_summary")
    p.add_argument("ultraplex_dir")
    p.add_argument("--prefix", default="upx")
    p.add_argument("--collapse-swapped", action="store_true")
    args = p.parse_args()

    total, assigned, unassigned, sm = load_seqmux(Path(args.seqmux_summary))
    up = load_ultraplex(Path(args.ultraplex_dir), args.prefix)
    if args.collapse_swapped:
        up = collapse_swapped(up)

    up_sum = sum(up.values())
    print(f"seqmux total={total} assigned={assigned} unassigned={unassigned}")
    print(f"ultraplex assigned_files={len(up)} assigned_sum={up_sum}")

    errors: list[str] = []
    if assigned != up_sum:
        errors.append(f"assigned {assigned} != ultraplex sum {up_sum} (delta {assigned - up_sum})")

    all_names = sorted(set(sm) | set(up))
    for name in all_names:
        a = sm.get(name, 0)
        b = up.get(name, 0)
        if a != b:
            errors.append(f"{name}: seqmux={a} ultraplex={b} delta={a - b}")

    if errors:
        print(f"DIFF {len(errors)}:")
        for e in errors[:80]:
            print(" ", e)
        return 1
    print("OK per-sample counts match")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
