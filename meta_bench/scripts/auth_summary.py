#!/usr/bin/env python3
"""Ground-truth authentication summary for the aMeta metagenomic benchmark.

This is a *post-hoc* analysis tool, deliberately kept out of the Snakemake
pipeline because it needs the simulation ground truth (which the pipeline never
sees). Given the aMeta AUTHENTICATION outputs and the ground-truth composition
tables, it produces, for every profiler:

  1. A heatmap with the same layout as aMeta's overview_heatmap_scores (taxa x
     samples, cell text = authentication score) but recoloured by ground truth:
       - taxon NOT in the sample  -> red   (i.e. a false positive)
       - MODERN taxon in sample   -> blue
       - ANCIENT taxon in sample  -> green
     Colour intensity encodes the taxon's abundance in that sample.

  2. An overall score difference vs the "perfect" score (0 for taxa not in the
     sample, --ancient-score for ancient genomes, --modern-score for modern).
     Two numbers per profiler:
       - full           : every ground-truth taxon counts.
       - kraken_adjusted : ground-truth taxa that KrakenUniq never surfaced for
                           authentication are dropped, so a KrakenUniq miss does
                           not penalise the downstream profiler.

  3. Detection metrics per sample (and their mean across samples):
       - jaccard_all         : Jaccard of {taxa found with score > 0} vs {taxa in
                               the sample}.
       - jaccard_ancient     : same, restricted to ancient genomes.
       - mean_reads_per_node : mean aligned reads (TotalAlignmentsOnReference)
                               over the sample's authenticated nodes.
       - evenness_tp/evenness_fp : mean breadth of coverage (fraction of the
                               reference covered) over found nodes that are true
                               positives / false positives.

Taxa are matched between the profiler output, the ground truth and KrakenUniq by
normalised genus+species name (the reference taxids in accessions.tsv are often
strain-level while the profilers report species-level, so a raw taxid join
misses several species).
"""

import argparse
import csv
import os
import re
import subprocess
import sys


# --------------------------------------------------------------------------- #
# Name handling
# --------------------------------------------------------------------------- #
def norm_key(name):
    """Normalise an organism/species label to a 'genus species' matching key."""
    if not name:
        return ""
    s = name.replace("_", " ").strip()
    if s.lower().endswith(".fa"):
        s = s[:-3]
    toks = re.findall(r"[A-Za-z0-9]+", s)
    if not toks:
        return ""
    return " ".join(toks[:2]).lower()


def display_from_basename(basename):
    """'Yersinia_pestis_CO92.fa' -> 'Yersinia pestis'."""
    s = basename[:-3] if basename.lower().endswith(".fa") else basename
    toks = re.findall(r"[A-Za-z0-9]+", s.replace("_", " "))
    if len(toks) >= 2:
        return toks[0].capitalize() + " " + toks[1].lower()
    return toks[0].capitalize() if toks else basename


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def load_accessions(path):
    """species_basename -> class ; also key -> class for name matching."""
    cls_by_key = {}
    host_keys = set()
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 2 or f[0] == "species":
                continue
            species, klass = f[0], f[1]
            key = norm_key(species)
            cls_by_key[key] = klass
            if klass == "host":
                host_keys.add(key)
    return cls_by_key, host_keys


def load_truth(path, category, samples):
    """Parse a ground-truth abundance table.

    Rows: <species>.fa\tab_sample1\tab_sample2...  columns are positional
    sample1..sampleN. Returns key -> {category, display, abund:{sample:frac}}.
    """
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            basename = f[0]
            key = norm_key(basename)
            abund = {}
            for i, val in enumerate(f[1:]):
                if i >= len(samples):
                    break
                try:
                    v = float(val)
                except ValueError:
                    v = 0.0
                if v > 0:
                    abund[samples[i]] = v
            if abund:
                out[key] = {
                    "category": category,
                    "display": display_from_basename(basename),
                    "abund": abund,
                }
    return out


def read_score_folder(folder):
    """Return (organism_display, score) for one AUTHENTICATION/<taxid> folder."""
    organism, score = "", 0
    sp = os.path.join(folder, "authentication_scores.txt")
    if os.path.isfile(sp):
        with open(sp) as fh:
            rows = [r.rstrip("\n").split("\t") for r in fh if r.strip()]
        if len(rows) >= 2 and len(rows[1]) >= 2:
            organism = rows[1][0].strip()
            try:
                score = int(float(rows[1][1]))
            except ValueError:
                print("v error")
                score = 0
    if not organism:
        nl = os.path.join(folder, "node_list.txt")
        if os.path.isfile(nl):
            with open(nl) as fh:
                first = fh.readline().strip()
            organism = first
    # print(organism, score)
    return organism, score


def read_node_reads(folder):
    """Aligned reads on the node's reference (TotalAlignmentsOnReference), or 0.

    Handles both MaltExtract's layout (header carries the column, no rowname) and
    the ngsLCA layout (header omits the leading taxid rowname the data row carries).
    """
    rd = os.path.join(folder, "MaltExtract_output", "default", "readDist")
    if not os.path.isdir(rd):
        return 0
    for fn in os.listdir(rd):
        if not fn.endswith("_alignmentDist.txt"):
            continue
        with open(os.path.join(rd, fn)) as fh:
            rows = [r.rstrip("\n").split("\t") for r in fh if r.strip()]
        if len(rows) < 2:
            return 0
        header, data = rows[0], rows[1]
        vals = data[1:] if len(data) == len(header) + 1 else data
        cell = dict(zip(header, vals)).get("TotalAlignmentsOnReference", 0)
        try:
            return int(float(cell))
        except (ValueError, TypeError):
            return 0
    return 0


def node_breadth(folder):
    """Breadth of coverage: fraction of reference positions with depth > 0.

    This is the "X% of genome covered" figure aMeta prints on its evenness-of-
    coverage plot. Computed from the per-base `samtools depth -a` file
    (`breadth_of_coverage`, cols: ref, pos, depth) with awk, since those files are
    large (up to hundreds of MB). Returns a float in [0, 1], or None if unavailable.
    """
    path = os.path.join(folder, "breadth_of_coverage")
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None
    try:
        r = subprocess.run(
            ["awk", '{t++; if($3>0) c++} END{if(t>0) printf "%.6f", c/t}', path],
            capture_output=True, text=True)
    except OSError:
        return None
    out = r.stdout.strip()
    try:
        return float(out) if out else None
    except ValueError:
        return None


def load_profiler_scores(auth_dir, profiler, samples):
    """(sample, key) -> {display, score, taxid} for one profiler."""
    scores = {}
    base = os.path.join(auth_dir, profiler)
    for sample in samples:
        sdir = os.path.join(base, sample)
        if not os.path.isdir(sdir):
            continue
        for taxid in sorted(os.listdir(sdir)):
            folder = os.path.join(sdir, taxid)
            if not os.path.isdir(folder):
                continue
            organism, score = read_score_folder(folder)
            reads = read_node_reads(folder)
            key = norm_key(organism) or taxid
            display = organism or taxid
            prev = scores.get((sample, key))
            if prev is None or score > prev["score"]:
                scores[(sample, key)] = {
                    "display": display,
                    "score": score,
                    "taxid": taxid,
                    "reads": reads,
                }
    return scores


def load_kraken_detected(results_dir, samples):
    """sample -> set of normalised keys that KrakenUniq surfaced (taxID.species)."""
    detected = {s: set() for s in samples}
    kdir = os.path.join(results_dir, "KRAKENUNIQ")
    for sample in samples:
        sp_file = os.path.join(kdir, sample, "taxID.species")
        filt = os.path.join(kdir, sample, "krakenuniq.output.filtered")
        if not os.path.isfile(sp_file):
            continue
        taxids = set()
        with open(sp_file) as fh:
            for line in fh:
                t = line.strip()
                if t:
                    taxids.add(t)
        name_by_taxid = {}
        if os.path.isfile(filt):
            with open(filt) as fh:
                for line in fh:
                    if line.startswith("#") or line.startswith("%"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 9:
                        name_by_taxid[f[6].strip()] = f[8].strip()
        for t in taxids:
            nm = name_by_taxid.get(t, "")
            key = norm_key(nm) if nm else t
            detected[sample].add(key)
    return detected


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
CAT_COLOR = {
    "ancient": (0.15, 0.65, 0.20),   # green
    "modern": (0.20, 0.45, 0.85),    # blue
    "absent": (0.85, 0.20, 0.20),    # red
}


def blend(base, intensity):
    """White -> base colour by intensity in [0,1], floored so it stays visible."""
    a = 0.30 + 0.70 * max(0.0, min(1.0, intensity))
    return tuple(1.0 * (1 - a) + c * a for c in base)


def luminance(rgb):
    return 0.299 * rgb[0] + 0.587 * rgb[1] + 0.750 * rgb[2]


def build_matrix(profiler_scores, truth, samples, include_host, host_keys):
    """Assemble ordered rows and per-cell records for one profiler."""
    keys = set()
    for (_s, key) in profiler_scores:
        keys.add(key)
    for key, info in truth.items():
        if any(s in info["abund"] for s in samples):
            keys.add(key)
    if not include_host:
        keys -= host_keys

    # canonical display label per key: ground truth name wins, else profiler name
    label = {}
    for key in keys:
        if key in truth:
            label[key] = truth[key]["display"]
    for (_s, key), rec in profiler_scores.items():
        if key in keys and key not in label:
            label[key] = rec["display"]
    for key in keys:
        label.setdefault(key, key)

    rows = sorted(keys, key=lambda k: label[k].lower())

    # per-sample max abundance for intensity normalisation
    smax = {}
    for s in samples:
        vals = [info["abund"].get(s, 0.0) for info in truth.values()]
        smax[s] = max(vals) if vals and max(vals) > 0 else 1.0

    return rows, label, smax


def render_heatmap(profiler, rows, label, samples, profiler_scores, truth,
                   smax, detected, out_pdf, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    nrow, ncol = len(rows), len(samples)
    if nrow == 0:
        print(f"  [{profiler}] no taxa to plot, skipping heatmap")
        return

    fig_w = max(4.0, 2.2 + 1.1 * ncol)
    fig_h = max(3.0, 0.34 * nrow + 1.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for r, key in enumerate(rows):
        gt = truth.get(key)
        for c, sample in enumerate(samples):
            in_sample = gt is not None and sample in gt["abund"]
            if in_sample:
                category = gt["category"]
                intensity = gt["abund"][sample] / smax[sample]
                color = blend(CAT_COLOR[category], intensity)
                authed = (sample, key) in profiler_scores
                if authed:
                    txt = str(profiler_scores[(sample, key)]["score"])
                elif key not in detected.get(sample, set()):
                    # KrakenUniq never surfaced this taxon, so the profiler had
                    # no chance to authenticate it: leave the score blank but
                    # keep the ground-truth colour. Distinguishes a KrakenUniq
                    # recall miss from a profiler authentication miss (score 0).
                    txt = ""
                else:
                    # KrakenUniq surfaced it but authentication scored 0 here
                    # (a genuine profiler miss).
                    txt = "0"
            else:
                authed = (sample, key) in profiler_scores
                if authed:  # false positive
                    color = CAT_COLOR["absent"]
                    txt = str(profiler_scores[(sample, key)]["score"])
                else:       # true negative: leave blank
                    color = (1.0, 1.0, 1.0)
                    txt = ""
            ax.add_patch(plt.Rectangle((c, nrow - 1 - r), 1, 1,
                                       facecolor=color, edgecolor="white", lw=1.0))
            if txt != "":
                tc = "white" if luminance(color) < 0.55 else "black"
                ax.text(c + 0.5, nrow - 1 - r + 0.5, txt, ha="center",
                        va="center", fontsize=8, color=tc)

    ax.set_xlim(0, ncol)
    ax.set_ylim(0, nrow)
    ax.set_xticks([c + 0.5 for c in range(ncol)])
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    ax.set_yticks([nrow - 1 - r + 0.5 for r in range(nrow)])
    ax.set_yticklabels([label[k] for k in rows], fontsize=8)
    ax.set_xticks(range(ncol + 1), minor=True)
    ax.set_yticks(range(nrow + 1), minor=True)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(f"Authentication vs ground truth — {profiler}", fontsize=11)

    legend = [
        Patch(facecolor=CAT_COLOR["ancient"], label="ancient (in sample)"),
        Patch(facecolor=CAT_COLOR["modern"], label="modern (in sample)"),
        Patch(facecolor=CAT_COLOR["absent"], label="not in sample (false positive)"),
        Patch(facecolor=(0.82, 0.82, 0.82), label="intensity ~ abundance"),
        Patch(facecolor=(0.82, 0.82, 0.82), edgecolor="0.4", hatch="////",
              label="coloured + blank = KrakenUniq miss"),
    ]
    # place the legend a fixed ~0.7in below the axis so it clears the rotated
    # sample tick labels regardless of figure height
    offset = -0.7 / fig_h
    ax.legend(handles=legend, bbox_to_anchor=(0.5, offset), loc="upper center",
              ncol=2, fontsize=8, frameon=False)

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [{profiler}] heatmap -> {out_pdf}")


# --------------------------------------------------------------------------- #
# Score difference
# --------------------------------------------------------------------------- #
def score_differences(samples, profiler_scores, truth, detected,
                      ancient_score, modern_score, include_host, host_keys):
    """Return per-sample (full, kraken_adjusted) L1 differences vs perfect."""
    perfect = {"ancient": ancient_score, "modern": modern_score}
    per_sample = {}
    for sample in samples:
        full = 0
        adj = 0
        # union of authenticated taxa and in-sample ground-truth taxa
        keys = {k for (s, k) in profiler_scores if s == sample}
        for key, info in truth.items():
            if sample in info["abund"]:
                keys.add(key)
        if not include_host:
            keys -= host_keys
        for key in keys:
            gt = truth.get(key)
            in_sample = gt is not None and sample in gt["abund"]
            want = perfect[gt["category"]] if in_sample else 0
            got = profiler_scores[(sample, key)]["score"] if (sample, key) in profiler_scores else 0
            d = abs(want - got)
            full += d
            # kraken_adjusted: drop in-sample taxa KrakenUniq never surfaced
            if in_sample and key not in detected.get(sample, set()):
                continue
            adj += d
        per_sample[sample] = (full, adj)
    return per_sample


# --------------------------------------------------------------------------- #
# Detection metrics
# --------------------------------------------------------------------------- #
def _jaccard(a, b):
    """|a n b| / |a u b|; NaN when both sets are empty (undefined)."""
    union = a | b
    return (len(a & b) / len(union)) if union else float("nan")


def _nanmean(vals):
    xs = [v for v in vals if v == v]  # drop NaN
    return (sum(xs) / len(xs)) if xs else float("nan")


def profiler_metrics(samples, profiler_scores, truth, include_host, host_keys,
                     auth_dir, profiler):
    """Per-sample detection metrics for one profiler.

    Returns sample -> (jaccard_all, jaccard_ancient, mean_reads_per_node,
                       evenness_tp, evenness_fp):
      jaccard_all     : Jaccard of {taxa found with score > 0} vs {taxa in sample}.
      jaccard_ancient : same, restricted to ancient genomes (a found taxon counts
                        only if ground truth classes it ancient).
      mean_reads_per_node : mean TotalAlignmentsOnReference over the sample's
                        authenticated nodes.
      evenness_tp/fp  : mean breadth of coverage (fraction of the reference covered)
                        over the found nodes that are true positives (in sample) /
                        false positives (not in sample).
    Host taxa are excluded unless include_host.
    """
    per_sample = {}
    for sample in samples:
        found = {k for (s, k), rec in profiler_scores.items()
                 if s == sample and rec["score"] > 0}
        present = {k for k, info in truth.items() if sample in info["abund"]}
        if not include_host:
            found -= host_keys
            present -= host_keys
        anc_found = {k for k in found
                     if k in truth and truth[k]["category"] == "ancient"}
        anc_present = {k for k in present if truth[k]["category"] == "ancient"}
        reads = [rec["reads"] for (s, k), rec in profiler_scores.items()
                 if s == sample and (include_host or k not in host_keys)]
        # breadth of coverage over found nodes, split into TP (in sample) / FP
        tp_b, fp_b = [], []
        for k in found:
            folder = os.path.join(auth_dir, profiler, sample,
                                  profiler_scores[(sample, k)]["taxid"])
            b = node_breadth(folder)
            if b is None:
                continue
            (tp_b if k in present else fp_b).append(b)
        per_sample[sample] = (
            _jaccard(found, present),
            _jaccard(anc_found, anc_present),
            (sum(reads) / len(reads)) if reads else float("nan"),
            (sum(tp_b) / len(tp_b)) if tp_b else float("nan"),
            (sum(fp_b) / len(fp_b)) if fp_b else float("nan"),
        )
    return per_sample


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True,
                    help="aMeta results dir (contains AUTHENTICATION/, KRAKENUNIQ/)")
    ap.add_argument("--accessions", required=True, help="config/accessions.tsv")
    ap.add_argument("--ancient", required=True,
                    help="ground_truth_ancient_microbes.txt")
    ap.add_argument("--modern", required=True,
                    help="ground_truth_modern_microbes.txt")
    ap.add_argument("--outdir", required=True, help="output directory")
    ap.add_argument("--profilers", default="",
                    help="comma-separated profilers (default: all under AUTHENTICATION/)")
    ap.add_argument("--ancient-score", type=float, default=10,
                    help="perfect score for an ancient genome (default 10)")
    ap.add_argument("--modern-score", type=float, default=6,
                    help="perfect score for a modern genome (default 6, the "
                         "empirical ceiling; use 5 for the on-merit value that "
                         "excludes the spuriously-earned ancient edit-distance point)")
    ap.add_argument("--include-host", action="store_true",
                    help="keep host taxa (e.g. Homo sapiens); default drops them")
    args = ap.parse_args()

    auth_dir = os.path.join(args.results, "AUTHENTICATION")
    if not os.path.isdir(auth_dir):
        sys.exit(f"No AUTHENTICATION dir under {args.results}")

    if args.profilers:
        profilers = [p for p in args.profilers.split(",") if p]
    else:
        profilers = sorted(d for d in os.listdir(auth_dir)
                           if os.path.isdir(os.path.join(auth_dir, d)))

    # sample universe = union of sample dirs seen across profilers, natural order
    sample_set = set()
    for p in profilers:
        pdir = os.path.join(auth_dir, p)
        for s in os.listdir(pdir):
            if os.path.isdir(os.path.join(pdir, s)):
                sample_set.add(s)

    def sample_num(s):
        m = re.search(r"(\d+)$", s)
        return (int(m.group(1)) if m else 1 << 30, s)
    samples = sorted(sample_set, key=sample_num)
    # ground-truth table columns are positional sample1..sampleN
    truth_samples = [f"sample{i + 1}" for i in range(max(
        (sample_num(s)[0] for s in samples if sample_num(s)[0] < (1 << 30)),
        default=0))]
    truth_samples = truth_samples or samples

    _cls_by_key, host_keys = load_accessions(args.accessions)
    truth = {}
    truth.update(load_truth(args.ancient, "ancient", truth_samples))
    truth.update(load_truth(args.modern, "modern", truth_samples))
    detected = load_kraken_detected(args.results, samples)

    os.makedirs(args.outdir, exist_ok=True)
    metric_rows = []
    for profiler in profilers:
        pscores = load_profiler_scores(auth_dir, profiler, samples)
        rows, label, smax = build_matrix(pscores, truth, samples,
                                         args.include_host, host_keys)
        render_heatmap(
            profiler, rows, label, samples, pscores, truth, smax, detected,
            os.path.join(args.outdir, f"auth_summary_heatmap.{profiler}.pdf"),
            os.path.join(args.outdir, f"auth_summary_heatmap.{profiler}.png"))

        diffs = score_differences(samples, pscores, truth, detected,
                                  args.ancient_score, args.modern_score,
                                  args.include_host, host_keys)
        metr = profiler_metrics(samples, pscores, truth,
                                args.include_host, host_keys, auth_dir, profiler)
        for sample in samples:
            full, adj = diffs[sample]
            j_all, j_anc, mreads, etp, efp = metr[sample]
            metric_rows.append((profiler, sample, full, adj, j_all, j_anc,
                                mreads, etp, efp))
        # MEAN across samples for every column
        metric_rows.append((
            profiler, "MEAN",
            _nanmean([diffs[s][0] for s in samples]),
            _nanmean([diffs[s][1] for s in samples]),
            _nanmean([metr[s][0] for s in samples]),
            _nanmean([metr[s][1] for s in samples]),
            _nanmean([metr[s][2] for s in samples]),
            _nanmean([metr[s][3] for s in samples]),
            _nanmean([metr[s][4] for s in samples]),
        ))

    def fmt(x, nd):
        if isinstance(x, float) and x != x:  # NaN
            return "NA"
        return f"{x:.{nd}f}" if isinstance(x, float) else str(x)

    cols = ["profiler", "sample", "diff_full", "diff_kraken_adjusted",
            "jaccard_all", "jaccard_ancient", "mean_reads_per_node",
            "evenness_tp", "evenness_fp"]

    def render_row(r):
        prof, sample, full, adj, ja, jaa, mr, etp, efp = r
        return [prof, sample, fmt(full, 1), fmt(adj, 1),
                fmt(ja, 3), fmt(jaa, 3), fmt(mr, 1), fmt(etp, 3), fmt(efp, 3)]

    out_metrics = os.path.join(args.outdir, "auth_metrics.tsv")
    with open(out_metrics, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in metric_rows:
            w.writerow(render_row(r))
    print(f"\nMetrics -> {out_metrics}")
    print("\n" + "\t".join(cols))
    for r in metric_rows:
        print("\t".join(render_row(r)))


if __name__ == "__main__":
    main()
