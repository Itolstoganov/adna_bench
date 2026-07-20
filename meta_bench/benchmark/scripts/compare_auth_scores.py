#!/usr/bin/env python3
"""Compare aMeta authentication-score tables from two profiler runs.

aMeta writes results/overview_heatmap_scores.txt as an organism x sample matrix of
integer authentication scores (rows = ORGANISM, columns = sample). This reads the
table from two runs (e.g. aMeta-malt and aMeta-ngslca) and emits a long side-by-side
comparison plus a short summary, so the MALT and strobealign+ngsLCA results can be
looked at together.

Usage:
  compare_auth_scores.py --malt PATH --ngslca PATH [--out OUT.tsv] [--threshold N]

PATH may be a run directory (results/overview_heatmap_scores.txt is appended) or the
table file itself.
"""
import argparse
import sys
from pathlib import Path


def resolve(path: Path) -> Path:
    if path.is_dir():
        return path / "results" / "overview_heatmap_scores.txt"
    return path


def load(path: Path):
    """Return {(organism, sample): score}."""
    p = resolve(path)
    if not p.exists():
        sys.exit(f"missing score table: {p}")
    lines = p.read_text().splitlines()
    if not lines:
        return {}, []
    samples = lines[0].split("\t")
    scores = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        organism = fields[0]
        for i, sample in enumerate(samples):
            if i + 1 >= len(fields):
                break
            val = fields[i + 1].strip()
            try:
                scores[(organism, sample)] = int(float(val))
            except ValueError:
                continue
    return scores, samples


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--malt", required=True, type=Path)
    ap.add_argument("--ngslca", required=True, type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--threshold", type=int, default=1,
                    help="score >= threshold counts as an authentic detection")
    args = ap.parse_args()

    malt, _ = load(args.malt)
    ngslca, _ = load(args.ngslca)

    keys = sorted(set(malt) | set(ngslca))
    rows = ["organism\tsample\tmalt_score\tngslca_score"]
    n_malt = n_ngslca = n_both = 0
    for organism, sample in keys:
        m = malt.get((organism, sample), 0)
        g = ngslca.get((organism, sample), 0)
        rows.append(f"{organism}\t{sample}\t{m}\t{g}")
        dm, dg = m >= args.threshold, g >= args.threshold
        n_malt += dm
        n_ngslca += dg
        n_both += dm and dg

    text = "\n".join(rows) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
        print(f"Wrote comparison ({len(keys)} organism x sample cells) to {args.out}",
              file=sys.stderr)
    else:
        sys.stdout.write(text)

    print(
        f"\n[summary] threshold>={args.threshold}: "
        f"MALT authentic={n_malt}, ngsLCA authentic={n_ngslca}, "
        f"both={n_both}, MALT-only={n_malt - n_both}, ngsLCA-only={n_ngslca - n_both}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
