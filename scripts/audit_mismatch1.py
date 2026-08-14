#!/usr/bin/env python3
"""Classify I464 pairs under exact vs mismatch=1 SeqMux rules."""
import csv
import gzip
from collections import Counter

BC_LEN = 8


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def load_samples(path):
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append(
                (
                    row["SampleNumber"],
                    row["Barcode1"].upper(),
                    row["Barcode2"].upper(),
                )
            )
    return rows


def best_match(r1, r2, samples, mm1, mm2):
    """Return (sample|None, dist, swapped|None, n_winners)."""
    best = None
    winners = []
    for name, b1, b2 in samples:
        if len(r1) < len(b1) or len(r2) < len(b2):
            continue
        d1 = hamming(r1[: len(b1)], b1)
        d2 = hamming(r2[: len(b2)], b2)
        if d1 <= mm1 and d2 <= mm2:
            d = d1 + d2
            winners.append((d, False, name, d1, d2))
        d1s = hamming(r2[: len(b1)], b1)
        d2s = hamming(r1[: len(b2)], b2)
        if d1s <= mm1 and d2s <= mm2:
            d = d1s + d2s
            winners.append((d, True, name, d1s, d2s))
    if not winners:
        return None, None, None, 0
    mind = min(w[0] for w in winners)
    at = [w for w in winners if w[0] == mind]
    names = {w[2] for w in at}
    if len(names) != 1:
        return None, mind, None, len(names)
    name = names.pop()
    swapped = all(w[1] for w in at if w[2] == name)
    return name, mind, swapped, 1


def exact_lookup(samples):
    canon = {}
    swap = {}
    bc1 = {}
    bc2 = {}
    for name, b1, b2 in samples:
        canon[(b1, b2)] = name
        swap[(b2, b1)] = name
        bc1.setdefault(b1, set()).add(name)
        bc2.setdefault(b2, set()).add(name)
    return canon, swap, bc1, bc2


def iter_pairs(r1p, r2p, max_reads=0):
    open1 = gzip.open if r1p.endswith(".gz") else open
    open2 = gzip.open if r2p.endswith(".gz") else open
    n = 0
    with open1(r1p, "rt") as f1, open2(r2p, "rt") as f2:
        while True:
            t1 = f1.readline()
            if not t1:
                break
            s1 = f1.readline()
            f1.readline()
            f1.readline()
            t2 = f2.readline()
            s2 = f2.readline()
            f2.readline()
            f2.readline()
            n += 1
            yield s1.strip().upper()[:BC_LEN], s2.strip().upper()[:BC_LEN]
            if max_reads and n >= max_reads:
                break


def main():
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--table")
    p.add_argument("--r1")
    p.add_argument("--r2")
    p.add_argument("--max-reads", type=int, default=0)
    p.add_argument("--out", default="-")
    args = p.parse_args()
    samples = load_samples(args.table)
    canon, swap, bc1map, bc2map = exact_lookup(samples)
    c = Counter()
    extra_by_sample = Counter()
    mis_from_to = Counter()

    for r1, r2 in iter_pairs(args.r1, args.r2, args.max_reads):
        c["total"] += 1
        exact = canon.get((r1, r2))
        swapped = False
        if exact is None:
            exact = swap.get((r1, r2))
            swapped = exact is not None
        mm_name, dist, mm_swapped, nw = best_match(r1, r2, samples, 1, 1)

        chimeric = False
        # Exact barcode pieces from two different samples
        s1 = bc1map.get(r1, set()) | bc2map.get(r1, set())
        s2 = bc1map.get(r2, set()) | bc2map.get(r2, set())
        if s1 and s2 and s1.isdisjoint(s2):
            chimeric = True
            c["chimeric_exact_pieces"] += 1

        if exact:
            c["exact"] += 1
            if swapped:
                c["exact_swapped"] += 1
            else:
                c["exact_canonical"] += 1
            if mm_name and mm_name != exact:
                c["mm1_overrides_exact"] += 1
        else:
            c["exact_unmatched"] += 1
            if nw > 1:
                c["mm1_ambiguous"] += 1
            elif mm_name is None:
                c["mm1_still_unmatched"] += 1
            else:
                c["mm1_new_assigned"] += 1
                extra_by_sample[mm_name] += 1
                if dist == 1:
                    c["mm1_new_dist1"] += 1
                elif dist == 2:
                    c["mm1_new_dist2"] += 1
                if chimeric:
                    c["mm1_absorbed_chimeric"] += 1
                    mis_from_to["chimeric->" + mm_name] += 1
                else:
                    c["mm1_recovered_errors"] += 1

    out = open(args.out, "w") if args.out != "-" else __import__("sys").stdout
    for k, v in sorted(c.items()):
        out.write(f"{k}\t{v}\n")
    out.write("\n# extra assignments by sample\n")
    for k, v in extra_by_sample.most_common():
        out.write(f"extra:{k}\t{v}\n")
    if args.out != "-":
        out.close()
    print("wrote", args.out, "total", c["total"], "mm1_new", c.get("mm1_new_assigned", 0))


if __name__ == "__main__":
    main()
