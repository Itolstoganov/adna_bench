#!/usr/bin/env python3
"""Set up a gargammel working directory and run gargammel.pl.

This is the single, de-duplicated home of the gargammel invocation shared by
sim_bench (top-level workflow/) and meta_bench (meta_bench/workflow/). It builds
the ``{workdir}/{endo,cont,bact}`` layout gargammel expects, then runs
``gargammel.pl`` with a caller-supplied option string.

Bacterial component, two modes:
  * ``--bact-dir DIR``            symlink DIR straight in as ``bact`` (sim_bench;
                                  gargammel auto-composes from all fastas).
  * ``--bact-list FILE`` +        create ``bact/`` and symlink each species FASTA
    ``--bact-genome-dir DIR``     named in the list from DIR, writing FILE to
                                  ``bact/list`` (meta_bench; explicit composition).

The fragment-size model, read length and damage flags are passed verbatim via
``--gargammel-args`` so each caller controls them (e.g. ``-l 40`` vs
``-rl 125 --loc .. --scale .. -damage ..``).
"""

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path


def symlink_relative(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    rel = os.path.relpath(src.resolve(), dst.parent.resolve())
    dst.symlink_to(rel)


def setup_component(paths, workdir: Path, name: str) -> None:
    """Symlink one or more FASTAs as {workdir}/{name}/{name}.{i}.fa."""
    comp_dir = workdir / name
    comp_dir.mkdir(parents=True, exist_ok=True)
    for i, p in enumerate(paths, start=1):
        symlink_relative(Path(p), comp_dir / f"{name}.{i}.fa")


def setup_bact_from_list(bact_list: Path, genome_dir: Path, workdir: Path) -> None:
    """Symlink each species FASTA named in bact_list into {workdir}/bact and
    copy the list to bact/list."""
    bact_dir = workdir / "bact"
    bact_dir.mkdir(parents=True, exist_ok=True)
    listed = []
    with open(bact_list) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            fname = line.split('\t')[0].split()[0]
            listed.append(fname)
            src = genome_dir / fname
            if not src.exists():
                sys.exit(f"Error: genome {src} referenced in {bact_list} not found")
            symlink_relative(src, bact_dir / fname)
    # Write the list verbatim so gargammel sees the exact abundances.
    (bact_dir / "list").write_text(Path(bact_list).read_text())


def main() -> None:
    ap = argparse.ArgumentParser(description="Run gargammel with a prepared workdir")
    ap.add_argument("--workdir", type=Path, required=True)
    ap.add_argument("--out", required=True, help="gargammel -o output prefix")
    ap.add_argument("--gargammel", required=True, help="path to gargammel.pl")
    ap.add_argument("--num", type=int, required=True, help="gargammel -n fragments")
    ap.add_argument("--comp", required=True, help="gargammel --comp bact,cont,endo")
    ap.add_argument("--endo", action="append", default=[], help="endogenous FASTA (repeatable)")
    ap.add_argument("--cont", action="append", default=[], help="contaminant FASTA (repeatable)")
    ap.add_argument("--bact-dir", default=None, help="symlink this dir as bact/ (sim_bench)")
    ap.add_argument("--bact-list", default=None, help="composition list -> bact/list (meta)")
    ap.add_argument("--bact-genome-dir", default=None, help="dir holding listed species FASTAs")
    ap.add_argument("--gargammel-args", default="", help="extra gargammel args (frag/rl/damage)")
    args = ap.parse_args()

    workdir = args.workdir
    workdir.mkdir(parents=True, exist_ok=True)

    if args.endo:
        setup_component(args.endo, workdir, "endo")
    if args.cont:
        setup_component(args.cont, workdir, "cont")

    if args.bact_dir:
        symlink_relative(Path(args.bact_dir), workdir / "bact")
    elif args.bact_list:
        if not args.bact_genome_dir:
            sys.exit("Error: --bact-list requires --bact-genome-dir")
        setup_bact_from_list(Path(args.bact_list), Path(args.bact_genome_dir), workdir)
    else:
        sys.exit("Error: provide either --bact-dir or --bact-list")

    # print(args.comp, sum([float(x) for x in args.comp]))
    comp = ','.join([str(float(x)) for x in args.comp.split(",")])
    cmd = (
        ["perl", args.gargammel, "-n", str(args.num), "--comp", comp]
        + shlex.split(args.gargammel_args)
        + ["-o", args.out, str(workdir)]
    )
    print("Running:", " ".join(shlex.quote(c) for c in cmd), file=sys.stderr)
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
