#!/usr/bin/env python3
"""Build a gargammel bacterial composition `list` for one sample/pass.

Reads a ground-truth abundance table (rows = ``species.fa<TAB>ab1<TAB>ab2...``,
one column per sample, as in
meta_bench/scripts/aMeta/ground_truth_{ancient,modern}_microbes.txt) and writes
the ``species.fa<TAB>abundance`` list gargammel expects for the given sample.

All species rows are emitted (absent ones with 0 abundance), mirroring the
original gargammel_sim.sh which fed gargammel the full list per pass.
"""

import argparse
import sys
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", required=True, type=Path, help="abundance table")
    ap.add_argument("--sample", required=True, type=int,
                    help="1-based sample number (selects table column)")
    ap.add_argument("--out", required=True, type=Path, help="output list path")
    args = ap.parse_args()

    rows = []
    with open(args.table) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            name = fields[0]
            if args.sample >= len(fields):
                sys.exit(f"Table {args.table} has no column for sample {args.sample} "
                         f"(row '{name}' has {len(fields) - 1} value columns)")
            abundance = fields[args.sample]
            rows.append((name, abundance))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as out:
        for name, abundance in rows:
            out.write(f"{name}\t{abundance}\n")

    present = sum(1 for _, a in rows if float(a) > 0)
    print(f"Wrote {len(rows)} species ({present} present) for sample {args.sample} "
          f"to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
