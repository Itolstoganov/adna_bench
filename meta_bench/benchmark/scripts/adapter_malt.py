#!/usr/bin/env python3
"""Adapter: aMeta MALT abundance matrix -> normalized taxon-abundance table.

The normalized contract consumed by score.py is a long TSV:
    sample <TAB> taxid <TAB> taxon <TAB> count

aMeta writes results/MALT_ABUNDANCE_MATRIX_SAM/malt_abundance_matrix_sam.txt as a
wide matrix: a header line of sample names, then one row per taxon
``<taxon name>\t<count sample1>\t<count sample2> ...``. MALT reports species
binomials (spaces) without taxids, so taxid is emitted as "NA" and score.py
falls back to name matching. (A future enrichment could recover taxids from the
MALT SAM RNAME ``accession|tax|taxid`` fields.)
"""

import argparse
import sys
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True, type=Path,
                    help="malt_abundance_matrix_sam.txt")
    ap.add_argument("--out", required=True, type=Path,
                    help="normalized taxon_abundance.tsv")
    args = ap.parse_args()

    lines = args.matrix.read_text().splitlines()
    if not lines:
        sys.exit(f"Empty matrix: {args.matrix}")

    samples = lines[0].split("\t")
    # Some aMeta versions leave a leading empty cell for the taxon column; drop it.
    if samples and samples[0] == "":
        samples = samples[1:]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(args.out, "w") as out:
        out.write("sample\ttaxid\ttaxon\tcount\n")
        for line in lines[1:]:
            if not line.strip():
                continue
            fields = line.split("\t")
            taxon = fields[0]
            counts = fields[1:]
            for i, sample in enumerate(samples):
                if i >= len(counts):
                    break
                val = counts[i].strip()
                if not val or val in ("0", "NA"):
                    continue
                out.write(f"{sample}\tNA\t{taxon}\t{val}\n")
                n += 1
    print(f"Wrote {n} (sample, taxon) abundance rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
