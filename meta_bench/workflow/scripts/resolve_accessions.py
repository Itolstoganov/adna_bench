#!/usr/bin/env python3
"""Snapshot accessions.tsv: fill empty accession/taxid columns by resolving each
species' reference genome via the NCBI ``datasets`` CLI, so the exact assemblies
are pinned for reproducibility. Rewrites the file in place (preserving comments).

    python meta_bench/workflow/scripts/resolve_accessions.py \
        --accessions meta_bench/config/accessions.tsv

Requires ncbi-datasets-cli. Rows that already have an accession are left as-is.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path


# RefSeq reference/complete first, then any RefSeq, then GenBank complete, then any
# GenBank — so taxa with no RefSeq assembly (e.g. Vermamoeba vermiformis) still resolve.
RESOLVE_ATTEMPTS = (
    ["--assembly-source", "RefSeq", "--reference"],
    ["--assembly-source", "RefSeq", "--assembly-level", "complete"],
    ["--assembly-source", "RefSeq"],
    ["--assembly-source", "GenBank", "--assembly-level", "complete"],
    ["--assembly-source", "GenBank"],
)


def resolve(name: str):
    """Return (accession, taxid) for a taxon name, or (None, None)."""
    for extra in RESOLVE_ATTEMPTS:
        cmd = ["datasets", "summary", "genome", "taxon", name,
               "--as-json-lines"] + extra
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        except subprocess.CalledProcessError:
            continue
        for line in out.splitlines():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            acc = rec.get("accession")
            taxid = str(rec.get("organism", {}).get("tax_id", "") or "")
            if acc:
                return acc, taxid
    return None, None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accessions", required=True, type=Path)
    args = ap.parse_args()

    lines = args.accessions.read_text().splitlines()
    out_lines = []
    for line in lines:
        if not line or line.startswith("#") or line.startswith("species\t"):
            out_lines.append(line)
            continue
        fields = line.split("\t")
        fields += [""] * (4 - len(fields))
        species, cls, taxid, acc = fields[:4]
        if not acc:
            name = species.replace("_", " ")
            print(f"Resolving {species} ({name})…", file=sys.stderr)
            r_acc, r_taxid = resolve(name)
            if r_acc:
                acc = r_acc
                taxid = taxid or r_taxid
                print(f"  -> {acc} (taxid {taxid})", file=sys.stderr)
            else:
                print(f"  !! could not resolve {species}", file=sys.stderr)
        out_lines.append("\t".join([species, cls, taxid, acc]).rstrip("\t"))

    args.accessions.write_text("\n".join(out_lines) + "\n")
    print(f"Updated {args.accessions}", file=sys.stderr)


if __name__ == "__main__":
    main()
