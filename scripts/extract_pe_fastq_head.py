#!/usr/bin/env python3
"""Copy the first N paired FASTQ records (plain or gzip) to gzip outputs."""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def open_text(path: Path, mode: str):
    if path.name.endswith(".gz") or "b" in mode:
        return gzip.open(path, mode + "t", compresslevel=1)
    return path.open(mode)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("-i", required=True, help="R1 FASTQ")
    p.add_argument("-I", dest="i2", required=True, help="R2 FASTQ")
    p.add_argument("-o", dest="o1", required=True, help="R1 output .fq.gz")
    p.add_argument("-O", dest="o2", required=True, help="R2 output .fq.gz")
    p.add_argument("-n", "--pairs", type=int, default=200_000)
    args = p.parse_args()

    n_lines = args.pairs * 4
    Path(args.o1).parent.mkdir(parents=True, exist_ok=True)
    with (
        open_text(Path(args.i), "r") as r1,
        open_text(Path(args.i2), "r") as r2,
        gzip.open(args.o1, "wt", compresslevel=1) as w1,
        gzip.open(args.o2, "wt", compresslevel=1) as w2,
    ):
        for i in range(n_lines):
            l1 = r1.readline()
            l2 = r2.readline()
            if not l1 or not l2:
                raise SystemExit(f"EOF after {i // 4} complete pairs")
            w1.write(l1)
            w2.write(l2)
    print(f"wrote {args.pairs} pairs -> {args.o1} {args.o2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
