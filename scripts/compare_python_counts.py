#!/usr/bin/env python3
"""Compare SeqMux summary.tsv sample counts to a Python sample_read_counts.csv."""
import argparse
import sys


def load_seqmux(path):
    total = assigned = unassigned = None
    orient_c = orient_s = None
    samples = {}
    with open(path) as f:
        next(f)
        for line in f:
            k, v = line.rstrip("\n").split("\t")
            if k == "total_reads":
                total = int(v)
            elif k == "assigned":
                assigned = int(v)
            elif k == "unassigned":
                unassigned = int(v)
            elif k == "orientation_canonical":
                orient_c = int(v)
            elif k == "orientation_swapped":
                orient_s = int(v)
            elif k.startswith("sample:"):
                samples[k.split(":", 1)[1]] = int(v)
    return total, assigned, unassigned, orient_c, orient_s, samples


def load_python_counts(path):
    samples = {}
    with open(path) as f:
        header = f.readline()
        key = "read_pairs" if "read_pairs" in header else header.strip().split(",")[-1]
        for line in f:
            parts = line.rstrip("\n").split(",")
            samples[parts[0]] = int(parts[-1])
    return samples


def main():
    p = argparse.ArgumentParser()
    p.add_argument("seqmux_summary")
    p.add_argument("python_counts")
    p.add_argument("--expect-total", type=int)
    p.add_argument("--expect-assigned", type=int)
    p.add_argument("--expect-canonical", type=int)
    p.add_argument("--expect-swapped", type=int)
    args = p.parse_args()

    total, assigned, unassigned, oc, os_, sm = load_seqmux(args.seqmux_summary)
    py = load_python_counts(args.python_counts)
    errors = []

    if args.expect_total is not None and total != args.expect_total:
        errors.append(f"total {total} != {args.expect_total}")
    if args.expect_assigned is not None and assigned != args.expect_assigned:
        errors.append(f"assigned {assigned} != {args.expect_assigned}")
    if args.expect_canonical is not None and oc != args.expect_canonical:
        errors.append(f"canonical {oc} != {args.expect_canonical}")
    if args.expect_swapped is not None and os_ != args.expect_swapped:
        errors.append(f"swapped {os_} != {args.expect_swapped}")

    py_sum = sum(py.values())
    if assigned != py_sum:
        errors.append(f"assigned {assigned} != python sum {py_sum}")

    all_names = sorted(set(sm) | set(py))
    for name in all_names:
        a = sm.get(name, 0)
        b = py.get(name, 0)
        if a != b:
            errors.append(f"{name}: seqmux={a} python={b} delta={a - b}")

    print(f"total={total} assigned={assigned} unassigned={unassigned}")
    print(f"orientation canonical={oc} swapped={os_}")
    print(f"python_samples={len(py)} seqmux_samples={len(sm)}")
    if errors:
        print(f"FAIL {len(errors)} differences:")
        for e in errors[:50]:
            print(" ", e)
        sys.exit(1)
    print("OK: per-sample counts match")


if __name__ == "__main__":
    main()
