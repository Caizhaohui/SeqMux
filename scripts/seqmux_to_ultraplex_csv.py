#!/usr/bin/env python3
"""Convert a SeqMux sample table to Ultraplex barcode CSV.

Ultraplex PE 3' barcodes are scored on reverse-complemented R2, from that
sequence's 3' end — i.e. reverse complement of R2 5'. To match SeqMux
canonical dual-barcode PE (Barcode1 @ R1 5', Barcode2 @ R2 5'), write:

    Barcode1,RC(Barcode2):SampleName

See docs/ULTRAPLEX_BENCH.md.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def load_seqmux_table(path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise SystemExit(f"no header in {path}")
        fields = {h.strip().lower().replace(" ", "_").replace("-", "_"): h for h in reader.fieldnames}

        def col(*aliases: str) -> str:
            for a in aliases:
                if a in fields:
                    return fields[a]
            raise SystemExit(f"missing column {aliases} in {path}")

        c_name = col("samplenumber", "sample_number", "sample", "sample_name", "name")
        c_bc1 = col("barcode1", "barcode_1", "bc1")
        c_bc2 = col("barcode2", "barcode_2", "bc2")
        for rec in reader:
            name = rec[c_name].strip()
            bc1 = rec[c_bc1].strip().upper()
            bc2 = rec[c_bc2].strip().upper()
            if not name or not bc1 or not bc2:
                raise SystemExit(f"empty sample/barcode in {path}: {rec}")
            rows.append((name, bc1, bc2))
    return rows


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("seqmux_csv")
    p.add_argument("-o", "--output", required=True)
    p.add_argument(
        "--no-rc-barcode2",
        action="store_true",
        help="write Barcode2 as-is (wrong for SeqMux PE 5'/5'; diagnostic only)",
    )
    p.add_argument(
        "--include-swapped",
        action="store_true",
        help="also emit Barcode2,RC(Barcode1):{name}_swapped rows",
    )
    args = p.parse_args()

    rows = load_seqmux_table(Path(args.seqmux_csv))
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for name, bc1, bc2 in rows:
        three = bc2 if args.no_rc_barcode2 else reverse_complement(bc2)
        lines.append(f"{bc1},{three}:{name}")
        if args.include_swapped:
            three_s = bc1 if args.no_rc_barcode2 else reverse_complement(bc1)
            lines.append(f"{bc2},{three_s}:{name}_swapped")
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {len(lines)} rows -> {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
