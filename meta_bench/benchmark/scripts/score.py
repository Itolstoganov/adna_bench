#!/usr/bin/env python3
"""Score a profiler's normalized taxon-abundance table against the simulation's
per-read ground truth.

Ground truth comes from the simulation's per-read truth tables
(``{sim_results}/{sample}/truth.tsv``); only bacterial reads (class == bact) are
the detection targets. A taxon is "truly present" in a sample if it has
>= --min-truth-reads simulated bacterial reads, and "detected" if the profiler
reports >= --detect-threshold count.

Matching uses a normalized "genus species" name (lowercased, first two tokens),
which is present on both sides and robust to strain suffixes like
Yersinia_pestis_CO92 vs "Yersinia pestis". (The MALT abundance matrix carries no
taxids, so taxid keying would never align with the taxid-bearing truth; a future
taxid-carrying adapter like ngsLCA could switch to taxid matching.)

Outputs:
  {out}/per_sample_scores.tsv   sample TP FP FN precision recall f1
  {out}/summary.tsv             micro-averaged totals + overall precision/recall/f1
  {out}/calls.tsv               per (sample, taxon) truth vs detected, for inspection
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path


def norm_name(name: str) -> str:
    toks = name.replace("_", " ").split()
    return " ".join(toks[:2]).lower()


def match_key(taxid: str, name: str) -> str:
    # Name-based only: the taxon name is present on both sides, whereas taxids are
    # not (MALT matrix has none). Keying on name keeps truth and detections aligned.
    return f"name:{norm_name(name)}"


def load_truth(sim_results: Path, samples, min_reads: int):
    """sample -> {match_key -> (label_taxon, read_count)} for bacterial reads."""
    truth = {}
    for sample in samples:
        path = sim_results / sample / "truth.tsv"
        if not path.exists():
            sys.exit(f"Missing truth table: {path}")
        counts = defaultdict(int)
        labels = {}
        with open(path) as f:
            header = f.readline().rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(header)}
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if fields[idx["class"]] != "bact":
                    continue
                taxon = fields[idx["taxon"]]
                taxid = fields[idx["taxid"]]
                key = match_key(taxid, taxon)
                counts[key] += 1
                labels.setdefault(key, taxon)
        truth[sample] = {
            k: (labels[k], c) for k, c in counts.items() if c >= min_reads
        }
    return truth


def load_detections(path: Path, threshold: float):
    """sample -> {match_key -> (label_taxon, count)}."""
    det = defaultdict(dict)
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            fields = line.rstrip("\n").split("\t")
            sample = fields[idx["sample"]]
            taxid = fields[idx["taxid"]]
            taxon = fields[idx["taxon"]]
            try:
                count = float(fields[idx["count"]])
            except ValueError:
                continue
            if count < threshold:
                continue
            det[sample][match_key(taxid, taxon)] = (taxon, count)
    return det


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sim-results", required=True, type=Path)
    ap.add_argument("--samples", required=True, nargs="+")
    ap.add_argument("--detections", required=True, type=Path,
                    help="normalized taxon_abundance.tsv from an adapter")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-truth-reads", type=int, default=1)
    ap.add_argument("--detect-threshold", type=float, default=1)
    args = ap.parse_args()

    truth = load_truth(args.sim_results, args.samples, args.min_truth_reads)
    det = load_detections(args.detections, args.detect_threshold)

    args.out.mkdir(parents=True, exist_ok=True)
    per_sample = []
    tot_tp = tot_fp = tot_fn = 0
    with open(args.out / "calls.tsv", "w") as calls:
        calls.write("sample\tkey\ttaxon\ttrue_reads\tdetected_count\tcall\n")
        for sample in args.samples:
            t = truth.get(sample, {})
            d = det.get(sample, {})
            keys = set(t) | set(d)
            tp = fp = fn = 0
            for key in sorted(keys):
                in_t = key in t
                in_d = key in d
                taxon = (t[key][0] if in_t else d[key][0])
                treads = t[key][1] if in_t else 0
                dcount = d[key][1] if in_d else 0
                if in_t and in_d:
                    call = "TP"; tp += 1
                elif in_d:
                    call = "FP"; fp += 1
                else:
                    call = "FN"; fn += 1
                calls.write(f"{sample}\t{key}\t{taxon}\t{treads}\t{dcount}\t{call}\n")
            prec = tp / (tp + fp) if (tp + fp) else 0.0
            rec = tp / (tp + fn) if (tp + fn) else 0.0
            f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
            per_sample.append((sample, tp, fp, fn, prec, rec, f1))
            tot_tp += tp; tot_fp += fp; tot_fn += fn

    with open(args.out / "per_sample_scores.tsv", "w") as f:
        f.write("sample\tTP\tFP\tFN\tprecision\trecall\tf1\n")
        for s, tp, fp, fn, p, r, f1 in per_sample:
            f.write(f"{s}\t{tp}\t{fp}\t{fn}\t{p:.4f}\t{r:.4f}\t{f1:.4f}\n")

    mprec = tot_tp / (tot_tp + tot_fp) if (tot_tp + tot_fp) else 0.0
    mrec = tot_tp / (tot_tp + tot_fn) if (tot_tp + tot_fn) else 0.0
    mf1 = 2 * mprec * mrec / (mprec + mrec) if (mprec + mrec) else 0.0
    with open(args.out / "summary.tsv", "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"total_TP\t{tot_tp}\ntotal_FP\t{tot_fp}\ntotal_FN\t{tot_fn}\n")
        f.write(f"micro_precision\t{mprec:.4f}\n")
        f.write(f"micro_recall\t{mrec:.4f}\n")
        f.write(f"micro_f1\t{mf1:.4f}\n")

    print(f"TP={tot_tp} FP={tot_fp} FN={tot_fn}  "
          f"precision={mprec:.3f} recall={mrec:.3f} f1={mf1:.3f}", file=sys.stderr)
    print(f"Wrote scores to {args.out}/", file=sys.stderr)


if __name__ == "__main__":
    main()
