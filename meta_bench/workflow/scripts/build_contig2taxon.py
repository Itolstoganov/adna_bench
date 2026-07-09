#!/usr/bin/env python3
"""Build a contig -> (taxon, taxid) map from the downloaded genomes.

gargammel/fragSim encodes the source *contig* id as a bacterial read's ref_id.
Downloaded RefSeq genomes use accession contig ids (e.g. NC_003143.1), so to
label reads by species we map every contig header in every genome FASTA to the
species (the FASTA basename) and its taxid (from accessions.tsv).

Output TSV (consumed by shared/scripts/parse_gargammel.py --contig2taxon):
    contig_id <TAB> species <TAB> taxid
"""

import argparse
import sys
from pathlib import Path


def load_taxids(accessions: Path) -> dict:
    taxids = {}
    with open(accessions) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] == "species":  # header
                continue
            species = fields[0]
            taxid = fields[2] if len(fields) > 2 and fields[2] else "NA"
            taxids[species] = taxid
    return taxids


def contig_ids(fasta: Path):
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                yield line[1:].split()[0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genomes", nargs="+", required=True, type=Path,
                    help="genome FASTA files ({species}.fa)")
    ap.add_argument("--accessions", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    taxids = load_taxids(args.accessions)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(args.out, "w") as out:
        for fasta in args.genomes:
            species = fasta.name
            if species.endswith(".fa"):
                species = species[:-3]
            elif species.endswith(".fasta"):
                species = species[:-6]
            taxid = taxids.get(species, "NA")
            for contig in contig_ids(fasta):
                out.write(f"{contig}\t{species}\t{taxid}\n")
                n += 1
    print(f"Wrote {n} contig mappings to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
