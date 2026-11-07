#!/usr/bin/env python3
"""Download one reference genome for the metagenomic simulation.

Given a species (FASTA basename) and an optional pinned accession, fetch the
genome and write it to ``{out}`` (a single FASTA named ``{species}.fa``). When
no accession is pinned, resolve the reference / representative *complete* genome
by species name.

Uses the NCBI ``datasets`` CLI (bioconda: ncbi-datasets-cli) for assembly
accessions (GCF_/GCA_). Any other accession is taken as a single nucleotide
record (e.g. NC_012920.1, human chrM for the test config) and fetched from
NCBI E-utilities. The resolved accession is written to ``{out}.accession`` for
provenance / snapshotting.
"""

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path


def have(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def species_to_name(species: str) -> str:
    """`Yersinia_pestis_CO92` -> `Yersinia pestis` (drop strain suffix heuristically)."""
    return species.replace("_", " ")


# Preference order: a RefSeq reference/complete genome first, then any RefSeq,
# then a GenBank complete genome, then any GenBank assembly. The GenBank fallback
# matters for taxa with no RefSeq assembly (e.g. some protists like Vermamoeba
# vermiformis), which would otherwise resolve to nothing.
RESOLVE_ATTEMPTS = (
    ["--assembly-source", "RefSeq", "--reference"],
    ["--assembly-source", "RefSeq", "--assembly-level", "complete"],
    ["--assembly-source", "RefSeq"],
    ["--assembly-source", "GenBank", "--assembly-level", "complete"],
    ["--assembly-source", "GenBank"],
)


def resolve_accession(name: str) -> str:
    """Return the best available assembly accession for a taxon name."""
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
            if acc:
                return acc
    sys.exit(f"Could not resolve an accession for '{name}'. Pin one in accessions.tsv.")


def download(accession: str, out: Path) -> None:
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        zip_path = tdp / "genome.zip"
        subprocess.run(
            ["datasets", "download", "genome", "accession", accession,
             "--include", "genome", "--filename", str(zip_path)],
            check=True,
        )
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tdp)
        fnas = sorted((tdp / "ncbi_dataset" / "data").rglob("*.fna"))
        if not fnas:
            sys.exit(f"No FASTA found in download for {accession}")
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "wb") as w:
            for fna in fnas:
                with open(fna, "rb") as r:
                    shutil.copyfileobj(r, w)


EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
          "?db=nuccore&rettype=fasta&retmode=text&id={}")


def is_assembly(accession: str) -> bool:
    return accession.startswith(("GCF_", "GCA_"))


def download_nucleotide(accession: str, out: Path) -> None:
    with urllib.request.urlopen(EFETCH.format(accession), timeout=120) as r:
        fasta = r.read()
    if not fasta.startswith(b">"):
        sys.exit(f"efetch returned no FASTA for {accession}")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(fasta)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--species", required=True)
    ap.add_argument("--accession", default="", help="pinned accession ('' to resolve)")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    accession = args.accession.strip()
    if accession and not is_assembly(accession):
        download_nucleotide(accession, args.out)
        Path(str(args.out) + ".accession").write_text(accession + "\n")
        print(f"Wrote {args.out} ({accession})", file=sys.stderr)
        return

    if not have("datasets"):
        sys.exit("NCBI 'datasets' CLI not found. Install ncbi-datasets-cli "
                 "(bioconda) or pin accessions and use a different fetcher.")

    if not accession:
        name = species_to_name(args.species)
        print(f"Resolving reference genome for '{name}'…", file=sys.stderr)
        accession = resolve_accession(name)
        print(f"  -> {accession}", file=sys.stderr)

    download(accession, args.out)
    Path(str(args.out) + ".accession").write_text(accession + "\n")
    print(f"Wrote {args.out} ({accession})", file=sys.stderr)


if __name__ == "__main__":
    main()
