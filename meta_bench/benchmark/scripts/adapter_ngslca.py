#!/usr/bin/env python3
"""Adapter: aMeta ngsLCA abundance matrix -> normalized taxon-abundance table.

The normalized contract consumed by score.py is a long TSV:
    sample <TAB> taxid <TAB> taxon <TAB> count

The aligner+ngsLCA profiler writes
results/NGSLCA_ABUNDANCE_MATRIX/ngslca_abundance_matrix.txt as a wide matrix with
two leading columns then one count column per sample:
    taxid <TAB> taxon <TAB> <count sample1> <TAB> <count sample2> ...
Unlike MALT, ngsLCA reports NCBI taxids directly, so this adapter populates the
taxid column (score.py can then match on taxid rather than falling back to name).
"""

import argparse
import sys
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True, type=Path,
                    help="ngslca_abundance_matrix.txt")
    ap.add_argument("--out", required=True, type=Path,
                    help="normalized taxon_abundance.tsv")
    args = ap.parse_args()

    lines = args.matrix.read_text().splitlines()
    if not lines:
        sys.exit(f"Empty matrix: {args.matrix}")

    header = lines[0].split("\t")
    # header: taxid, taxon, <sample1>, <sample2>, ...
    samples = header[2:]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(args.out, "w") as out:
        out.write("sample\ttaxid\ttaxon\tcount\n")
        for line in lines[1:]:
            if not line.strip():
                continue
            fields = line.split("\t")
            taxid = fields[0]
            taxon = fields[1]
            counts = fields[2:]
            for i, sample in enumerate(samples):
                if i >= len(counts):
                    break
                val = counts[i].strip()
                if not val or val in ("0", "NA"):
                    continue
                out.write(f"{sample}\t{taxid}\t{taxon}\t{val}\n")
                n += 1
    print(f"Wrote {n} (sample, taxon) abundance rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
