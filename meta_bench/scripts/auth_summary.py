#!/usr/bin/env python3
"""Ground-truth authentication summary for the aMeta metagenomic benchmark.

A *post-hoc* analysis tool, deliberately kept out of the Snakemake pipeline
because it needs the simulation ground truth (which the pipeline never sees).
Given the aMeta AUTHENTICATION outputs and the per-read simulation truth, it
produces, for every profiler:

  1. A heatmap with the same layout as aMeta's overview_heatmap_scores (taxa x
     samples, cell text = authentication score) recoloured by ground truth:
     green = ancient in sample, blue = modern in sample, red = not in sample
     (a false positive), intensity ~ abundance. A cell that is white with a "0"
     was screened and correctly rejected; white and empty was never a candidate.
  2. `auth_metrics.tsv`: per sample plus a MEAN row -- `diff_full` and
     `diff_kraken_adjusted` (L1 distance from the perfect score, the latter
     ignoring taxa KrakenUniq never surfaced), `jaccard_all`/`jaccard_ancient`
     (found with score > 0 vs in-sample truth), `mean_reads_per_node`, and
     `evenness_tp`/`evenness_fp` (breadth of coverage over true/false positives).
  3. `auth_calls.tsv`: one row per (profiler, sample, taxon) with read support,
     score and a TP/FP/FN verdict -- the per-cell backing for both of the above.

Ground truth is the simulated *reads* (`results/sim/<sample>/truth.tsv`), not the
composition tables in meta_bench/config/, because gargammel pairs the `bact/list`
weights with the genome files in readdir order (gargammel.pl:942 vs :1340) and so
simulates a permutation of what those tables ask for; truth.tsv is derived from
the read names and is downstream of that bug. Its `class`/`pass` columns give the
category (bact+ancient -> ancient microbe, bact+modern -> modern, endo/cont ->
host). Truth, profiler and KrakenUniq are joined on the enclosing species taxid
via the NCBI taxonomy, which is required (`--taxonomy`, auto-detected under
`<results>/NGSLCA_DB`) and has no name-based fallback: a raw taxid join breaks
because 11 of this benchmark's 36 reference genomes are pinned below species rank
(Mycobacterium_avium = 439334, subsp. hominissuis, vs the profilers' 1764), and a
name join breaks because one taxid carries several binomials (43767 is
`Rhodococcus hoagii` in accessions.tsv, `Rhodococcus equi` in the DB) -- each
failure mode splitting one organism into a phantom false positive plus a phantom
false negative.
"""

import argparse
import csv
import os
import re
import subprocess
import sys

import pandas as pd

# gargammel's -damage deaminates the bacterial and endogenous pools but leaves
# present-day contamination (cont) undamaged, even in the ancient pass.
DAMAGED_CLASSES = ("bact", "endo")
HOST_CLASSES = ("endo", "cont")
TRUTH_COLS = ("class", "taxon", "taxid", "pass")

VERDICTS = ("TP", "FP", "FN_profiler", "FN_krakenuniq", "TN")

# (column, decimal places) for auth_metrics.tsv; the MEAN row and the formatting
# are both driven off this, so adding a metric here is the only edit needed.
METRIC_SPEC = (
    ("diff_full", 1),
    ("diff_kraken_adjusted", 1),
    ("jaccard_all", 3),
    ("jaccard_ancient", 3),
    ("mean_reads_per_node", 1),
    ("evenness_tp", 3),
    ("evenness_fp", 3),
)
METRIC_COLS = ["profiler", "sample"] + [name for name, _ in METRIC_SPEC]

CALL_COLS = ["profiler", "sample", "taxon", "category", "truth_reads",
             "truth_frac", "score", "aligned_reads", "krakenuniq_surfaced",
             "verdict"]


# --------------------------------------------------------------------------- #
# Taxonomy
# --------------------------------------------------------------------------- #
def dmp_rows(path):
    """Yield the stripped fields of each row of an NCBI `*.dmp` file."""
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="|"):
            yield [field.strip() for field in row]


class Taxonomy:
    """NCBI taxonomy, used to key every source on its enclosing species taxid."""

    def __init__(self, nodes_dmp, names_dmp=None, merged_dmp=None):
        self.parent, self.rank, self.sci, self.merged = {}, {}, {}, {}
        self._cache = {}
        self.unresolved = {}          # taxid -> label seen, for reporting
        for f in dmp_rows(nodes_dmp):
            if len(f) >= 3:
                self.parent[f[0]], self.rank[f[0]] = f[1], f[2]
        if names_dmp and os.path.isfile(names_dmp):
            for f in dmp_rows(names_dmp):
                if len(f) >= 4 and f[3] == "scientific name":
                    self.sci[f[0]] = f[1]
        # merged.dmp is optional and often absent from an aMeta DB dir; without
        # it a taxid retired between DB build and genome download stays
        # unresolved rather than silently joining to nothing.
        if merged_dmp and os.path.isfile(merged_dmp):
            for f in dmp_rows(merged_dmp):
                if len(f) >= 2:
                    self.merged[f[0]] = f[1]

    def species(self, taxid):
        """The enclosing species taxid, or None if there is no species ancestor."""
        taxid = str(taxid).strip()
        if taxid in self._cache:
            return self._cache[taxid]
        t = self.merged.get(taxid, taxid)
        seen = set()
        out = None
        while t in self.parent and t not in seen and t != "1":
            if self.rank.get(t) == "species":
                out = t
                break
            seen.add(t)
            t = self.parent[t]
        self._cache[taxid] = out
        return out

    def key(self, taxid, label=""):
        """Join key for a taxid: its species node, else the taxid itself.

        A node above species (KrakenUniq can report e.g. 1649845, rank "species
        group") cannot be rolled up, so it keeps its own taxid and simply will
        not match anything -- which is the honest outcome, not a silent merge.
        """
        taxid = str(taxid).strip()
        sp = self.species(taxid)
        if sp:
            return sp
        self.unresolved[taxid] = label or self.sci.get(taxid, "")
        return taxid

    def name(self, key):
        return self.sci.get(key, "")


def find_taxonomy(spec, results_dir):
    """Locate nodes.dmp from --taxonomy, or auto-detect it under the results dir."""
    if spec:
        candidates = [spec, os.path.join(spec, "nodes.dmp")]
    else:
        candidates = [os.path.join(results_dir, db, sub, "nodes.dmp")
                      for db in ("NGSLCA_DB", "MALT_DB", "KRAKENUNIQ_DB")
                      for sub in ("", "taxonomy")]
    for c in candidates:
        if os.path.isfile(c) and os.path.basename(c) == "nodes.dmp":
            return c
    sys.exit(
        "No nodes.dmp found -- taxa are joined by species taxid and there is no "
        "name-based fallback.\nLooked in:\n  " + "\n  ".join(candidates) +
        "\nPass --taxonomy <dir containing nodes.dmp (and names.dmp)>.")


def open_taxonomy(spec, results_dir):
    """Load the taxonomy next to the discovered nodes.dmp, reporting what is there."""
    nodes = find_taxonomy(spec, results_dir)
    dmp_dir = os.path.dirname(nodes)
    print(f"Taxonomy: {nodes}")
    tax = Taxonomy(nodes,
                   os.path.join(dmp_dir, "names.dmp"),
                   os.path.join(dmp_dir, "merged.dmp"))
    if not tax.sci:
        print("  (no names.dmp beside it -- labels fall back to source names)")
    if not tax.merged:
        print("  (no merged.dmp beside it -- retired taxids will not resolve)")
    return tax


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def display_from_basename(basename):
    """'Yersinia_pestis_CO92.fa' -> 'Yersinia pestis'."""
    s = basename[:-3] if basename.lower().endswith(".fa") else basename
    toks = re.findall(r"[A-Za-z0-9]+", s.replace("_", " "))
    if len(toks) >= 2:
        return toks[0].capitalize() + " " + toks[1].lower()
    return toks[0].capitalize() if toks else basename


def truth_path_for(spec, sample):
    """Resolve the truth.tsv for one sample from --truth (dir or {sample} pattern)."""
    if "{sample}" in spec:
        return spec.format(sample=sample)
    return os.path.join(spec, sample, "truth.tsv")


def load_truth(spec, samples, tax):
    """Build the ground truth from the per-read `truth.tsv` of each sample.

    A read is ancient when it comes from the damaged pass *and* belongs to a
    class gargammel deaminates (bact, endo -- not cont, which stays undamaged
    even in the ancient pass); a taxon is ancient if at least half its reads are.
    Reads are summed per species taxid, so the two host labels
    (Homo_sapiens_ancient / _modern, both 9606) collapse into one row.

    Returns (truth, host_keys) where truth is
      key -> {category, display, abund:{sample: frac}, reads:{sample: n}}
    and `frac` is the taxon's share of that sample's reads.
    """
    truth, host_keys = {}, set()
    for sample in samples:
        path = truth_path_for(spec, sample)
        if not os.path.isfile(path):
            sys.exit(f"No truth.tsv for {sample} at {path}")
        try:
            df = pd.read_csv(path, sep="\t", usecols=list(TRUTH_COLS), dtype=str)
        except ValueError as exc:
            sys.exit(f"{path}: expected columns {', '.join(TRUTH_COLS)} ({exc})")
        total = len(df)
        if not total:
            continue
        # resolve each distinct taxid once, not once per read
        keys = {r.taxid: tax.key(r.taxid, r.taxon)
                for r in df[["taxid", "taxon"]].drop_duplicates().itertuples()}
        df["key"] = df["taxid"].map(keys)
        df["is_ancient"] = ((df["pass"] == "ancient")
                            & df["class"].isin(DAMAGED_CLASSES))
        df["is_host"] = df["class"].isin(HOST_CLASSES)
        grouped = df.groupby("key").agg(reads=("taxid", "size"),
                                        ancient=("is_ancient", "sum"),
                                        host=("is_host", "any"),
                                        taxon=("taxon", "first"))
        for key, row in grouped.iterrows():
            if row.host:
                host_keys.add(key)
            info = truth.setdefault(key, {
                "category": "modern",
                "display": display_from_basename(row.taxon),
                "abund": {}, "reads": {},
            })
            info["reads"][sample] = int(row.reads)
            info["abund"][sample] = row.reads / total
            # category is a property of the genome, not the sample; ancient wins
            # as soon as any sample simulated it in the damaged pass.
            if row.ancient * 2 >= row.reads:
                info["category"] = "ancient"
    return truth, host_keys


def read_score_folder(folder):
    """Return (organism_display, score) for one AUTHENTICATION/<taxid> folder."""
    organism, score = "", 0
    path = os.path.join(folder, "authentication_scores.txt")
    if os.path.isfile(path):
        with open(path) as fh:
            rows = [r for r in csv.reader(fh, delimiter="\t") if any(r)]
        if len(rows) >= 2 and len(rows[1]) >= 2:
            organism = rows[1][0].strip()
            try:
                score = int(float(rows[1][1]))
            except ValueError:
                score = 0
    if not organism:
        # a node aMeta screened but never scored still names itself here
        node_list = os.path.join(folder, "node_list.txt")
        if os.path.isfile(node_list):
            with open(node_list) as fh:
                organism = fh.readline().strip()
    return organism, score


def read_node_reads(folder):
    """Aligned reads on the node's reference (TotalAlignmentsOnReference), or 0.

    Handles both MaltExtract's layout (header carries the column, no rowname) and
    the ngsLCA layout (header omits the leading taxid rowname the data row
    carries), so the columns cannot simply be zipped by a DictReader.
    """
    rd = os.path.join(folder, "MaltExtract_output", "default", "readDist")
    if not os.path.isdir(rd):
        return 0
    for fn in os.listdir(rd):
        if not fn.endswith("_alignmentDist.txt"):
            continue
        with open(os.path.join(rd, fn)) as fh:
            rows = [r for r in csv.reader(fh, delimiter="\t") if any(r)]
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


def load_profiler_scores(auth_dir, profiler, samples, tax):
    """(sample, key) -> {display, score, taxid, reads} for one profiler.

    The AUTHENTICATION/<profiler>/<sample>/<taxid>/ directory name *is* the
    KrakenUniq taxid, so the profiler and KrakenUniq sides already share ids by
    construction; rolling both up to species only has to bridge the truth side.
    """
    scores = {}
    for sample in samples:
        sdir = os.path.join(auth_dir, profiler, sample)
        if not os.path.isdir(sdir):
            continue
        for taxid in sorted(os.listdir(sdir)):
            folder = os.path.join(sdir, taxid)
            if not os.path.isdir(folder):
                continue
            organism, score = read_score_folder(folder)
            key = tax.key(taxid, organism)
            prev = scores.get((sample, key))
            if prev is None or score > prev["score"]:
                scores[(sample, key)] = {
                    "display": organism or tax.name(key) or taxid,
                    "score": score,
                    "taxid": taxid,
                    "reads": read_node_reads(folder),
                }
    return scores


def load_kraken_detected(results_dir, samples, tax):
    """sample -> set of species keys that KrakenUniq surfaced (taxID.species)."""
    detected = {s: set() for s in samples}
    for sample in samples:
        kdir = os.path.join(results_dir, "KRAKENUNIQ", sample)
        sp_file = os.path.join(kdir, "taxID.species")
        if not os.path.isfile(sp_file):
            continue
        with open(sp_file) as fh:
            taxids = {line.strip() for line in fh if line.strip()}
        # names are only for the unresolved-taxid report, so a missing report is fine
        names = {}
        report = os.path.join(kdir, "krakenuniq.output.filtered")
        if os.path.isfile(report):
            with open(report) as fh:
                for f in csv.reader(fh, delimiter="\t"):
                    if len(f) >= 9 and not f[0].startswith(("#", "%")):
                        names[f[6].strip()] = f[8].strip()
        detected[sample] = {tax.key(t, names.get(t, "")) for t in taxids}
    return detected


# --------------------------------------------------------------------------- #
# Shared views over (truth, profiler_scores)
# --------------------------------------------------------------------------- #
def is_present(truth, key, sample):
    """Was this taxon simulated into this sample?"""
    info = truth.get(key)
    return info is not None and sample in info["abund"]


def sample_keys(sample, profiler_scores, truth, include_host, host_keys):
    """Taxa worth a verdict in this sample: authenticated or simulated."""
    keys = {k for (s, k) in profiler_scores if s == sample}
    keys |= {k for k, info in truth.items() if sample in info["abund"]}
    return keys if include_host else keys - host_keys


def verdict_for(in_sample, score, surfaced):
    """TP / FP / TN, splitting misses by whether KrakenUniq ever offered the taxon."""
    if in_sample:
        if score > 0:
            return "TP"
        return "FN_profiler" if surfaced else "FN_krakenuniq"
    return "FP" if score > 0 else "TN"


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #
def _jaccard(a, b):
    """|a n b| / |a u b|; NaN when both sets are empty (undefined)."""
    union = a | b
    return (len(a & b) / len(union)) if union else float("nan")


def _mean(vals):
    """Mean of the non-NaN values; NaN when there are none."""
    xs = [v for v in vals if v == v]
    return (sum(xs) / len(xs)) if xs else float("nan")


def sample_metrics(sample, profiler_scores, truth, detected, perfect,
                   include_host, host_keys, auth_dir, profiler):
    """The METRIC_SPEC metrics for one (profiler, sample), as a name -> value dict."""
    keys = sample_keys(sample, profiler_scores, truth, include_host, host_keys)
    surfaced = detected.get(sample, set())

    diff_full = diff_adjusted = 0
    for key in keys:
        in_sample = is_present(truth, key, sample)
        rec = profiler_scores.get((sample, key))
        want = perfect[truth[key]["category"]] if in_sample else 0
        d = abs(want - (rec["score"] if rec else 0))
        diff_full += d
        # kraken_adjusted: an in-sample taxon KrakenUniq never surfaced was never
        # the profiler's to find, so it does not count against it
        if not (in_sample and key not in surfaced):
            diff_adjusted += d

    found = {k for k in keys
             if (sample, k) in profiler_scores
             and profiler_scores[(sample, k)]["score"] > 0}
    present = {k for k in keys if is_present(truth, k, sample)}
    ancient = {k for k in truth if truth[k]["category"] == "ancient"}

    # breadth of coverage over found nodes, split into TP (in sample) / FP
    breadth = {True: [], False: []}
    for key in found:
        b = node_breadth(os.path.join(auth_dir, profiler, sample,
                                      profiler_scores[(sample, key)]["taxid"]))
        if b is not None:
            breadth[key in present].append(b)

    reads = [rec["reads"] for (s, k), rec in profiler_scores.items()
             if s == sample and (include_host or k not in host_keys)]

    return {
        "diff_full": diff_full,
        "diff_kraken_adjusted": diff_adjusted,
        "jaccard_all": _jaccard(found, present),
        "jaccard_ancient": _jaccard(found & ancient, present & ancient),
        "mean_reads_per_node": _mean(reads),
        "evenness_tp": _mean(breadth[True]),
        "evenness_fp": _mean(breadth[False]),
    }


def metric_rows(profiler, samples, per_sample):
    """Per-sample rows plus a MEAN row, driven off METRIC_SPEC."""
    rows = [{"profiler": profiler, "sample": s, **per_sample[s]} for s in samples]
    rows.append({"profiler": profiler, "sample": "MEAN",
                 **{name: _mean([per_sample[s][name] for s in samples])
                    for name, _ in METRIC_SPEC}})
    return rows


def call_rows(profiler, samples, profiler_scores, truth, detected, label,
              include_host, host_keys):
    """One row per (sample, taxon) that is either simulated or authenticated."""
    rows = []
    for sample in samples:
        keys = sample_keys(sample, profiler_scores, truth, include_host, host_keys)
        for key in sorted(keys, key=lambda k: label.get(k, k).lower()):
            in_sample = is_present(truth, key, sample)
            rec = profiler_scores.get((sample, key))
            score = rec["score"] if rec else 0
            surfaced = key in detected.get(sample, set())
            gt = truth[key] if in_sample else None
            rows.append({
                "profiler": profiler,
                "sample": sample,
                "taxon": label.get(key, key),
                "category": gt["category"] if gt else "absent",
                "truth_reads": gt["reads"].get(sample, 0) if gt else 0,
                "truth_frac": f"{gt['abund'][sample]:.5f}" if gt else "0",
                "score": score,
                "aligned_reads": rec["reads"] if rec else 0,
                "krakenuniq_surfaced": "yes" if surfaced else "no",
                "verdict": verdict_for(in_sample, score, surfaced),
            })
    return rows


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
CAT_COLOR = {
    "ancient": (0.15, 0.65, 0.20),   # green
    "modern": (0.20, 0.45, 0.85),    # blue
    "absent": (0.85, 0.20, 0.20),    # red
}
WHITE = (1.0, 1.0, 1.0)


def blend(base, intensity):
    """White -> base colour by intensity in [0,1], floored so it stays visible."""
    a = 0.30 + 0.70 * max(0.0, min(1.0, intensity))
    return tuple(1.0 * (1 - a) + c * a for c in base)


def luminance(rgb):
    return 0.299 * rgb[0] + 0.587 * rgb[1] + 0.750 * rgb[2]


def build_matrix(profiler_scores, truth, samples, include_host, host_keys, tax):
    """Assemble the ordered rows, their labels, and the intensity normalisation."""
    # every authenticated node earns a row, including one that scored 0 in every
    # sample -- those cells render as a white "0" (a correct rejection), which is
    # information, not an empty row.
    keys = {k for (_s, k) in profiler_scores}
    keys |= {k for k, info in truth.items()
             if any(s in info["abund"] for s in samples)}
    if not include_host:
        keys -= host_keys

    # canonical display label per key: ground truth name wins, else profiler name
    label = {k: truth[k]["display"] for k in keys if k in truth}
    for (_s, key), rec in profiler_scores.items():
        if key in keys:
            label.setdefault(key, rec["display"])
    for key in keys:
        label.setdefault(key, tax.name(key) or key)

    # per-sample max abundance for intensity normalisation, over plotted rows
    # only -- the host dwarfs every microbe, so including it when it is not
    # drawn would flatten the whole colour range.
    smax = {}
    for s in samples:
        vals = [info["abund"].get(s, 0.0) for k, info in truth.items() if k in keys]
        smax[s] = max(vals, default=0.0) or 1.0

    return sorted(keys, key=lambda k: label[k].lower()), label, smax


def cell_style(key, sample, truth, profiler_scores, detected, smax):
    """(facecolour, cell text) for one heatmap cell."""
    rec = profiler_scores.get((sample, key))
    if is_present(truth, key, sample):
        gt = truth[key]
        color = blend(CAT_COLOR[gt["category"]], gt["abund"][sample] / smax[sample])
        if rec:
            return color, str(rec["score"])
        # KrakenUniq never surfaced it, so the profiler had no chance to
        # authenticate it: keep the ground-truth colour but leave the score
        # blank, which distinguishes that from a profiler miss (an explicit 0).
        if key not in detected.get(sample, set()):
            return color, ""
        return color, "0"
    if rec and rec["score"] > 0:
        return CAT_COLOR["absent"], str(rec["score"])
    # Both remaining cases are true negatives, so no red -- but a node that was
    # screened and scored 0 is a correct rejection the profiler earned, whereas a
    # blank cell was never a candidate at all.
    return WHITE, "0" if rec else ""


def row_label(key, label):
    """'Mycobacterium avium (1764)' -- the key is the species taxid it joined on."""
    name = label.get(key, key)
    return name if name == key else f"{name} ({key})"


def render_heatmap(profiler, rows, label, samples, profiler_scores, truth,
                   smax, detected, out_pdf, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, Rectangle

    nrow, ncol = len(rows), len(samples)
    if nrow == 0:
        print(f"  [{profiler}] no taxa to plot, skipping heatmap")
        return

    fig_w = max(4.0, 2.2 + 1.1 * ncol)
    fig_h = max(3.0, 0.34 * nrow + 1.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for r, key in enumerate(rows):
        for c, sample in enumerate(samples):
            color, txt = cell_style(key, sample, truth, profiler_scores,
                                    detected, smax)
            ax.add_patch(Rectangle((c, nrow - 1 - r), 1, 1,
                                   facecolor=color, edgecolor="white", lw=1.0))
            if txt:
                ax.text(c + 0.5, nrow - 1 - r + 0.5, txt, ha="center", va="center",
                        fontsize=8,
                        color="white" if luminance(color) < 0.55 else "black")

    ax.set_xlim(0, ncol)
    ax.set_ylim(0, nrow)
    ax.set_xticks([c + 0.5 for c in range(ncol)])
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    ax.set_yticks([nrow - 1 - r + 0.5 for r in range(nrow)])
    ax.set_yticklabels([row_label(k, label) for k in rows], fontsize=8)
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
        Patch(facecolor=WHITE, edgecolor="0.6",
              label='white "0" = screened, correctly rejected'),
        Patch(facecolor=(0.82, 0.82, 0.82), label="intensity ~ abundance"),
        Patch(facecolor=(0.82, 0.82, 0.82), edgecolor="0.4", hatch="////",
              label="coloured + blank = KrakenUniq miss"),
    ]
    # place the legend a fixed ~0.7in below the axis so it clears the rotated
    # sample tick labels regardless of figure height
    ax.legend(handles=legend, bbox_to_anchor=(0.5, -0.7 / fig_h),
              loc="upper center", ncol=2, fontsize=8, frameon=False)

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [{profiler}] heatmap -> {out_pdf}")


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def fmt_cell(value, places):
    """Render one table cell: NaN as NA, floats to `places`, everything as-is."""
    if isinstance(value, float):
        return "NA" if value != value else f"{value:.{places}f}"
    return str(value)


def render_row(row, cols, places):
    return [fmt_cell(row.get(c, ""), places.get(c, 3)) for c in cols]


def write_tsv(path, cols, rows, places=None):
    """Write dict rows as a TSV (LF-terminated, so awk/cut behave)."""
    places = places or {}
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for row in rows:
            w.writerow(render_row(row, cols, places))


def print_tsv(cols, rows, places=None):
    places = places or {}
    print("\n" + "\t".join(cols))
    for row in rows:
        print("\t".join(render_row(row, cols, places)))


def report_unresolved(tax):
    if not tax.unresolved:
        return
    print(f"\n{len(tax.unresolved)} taxid(s) with no species ancestor "
          f"(kept as their own key, so they match nothing):")
    for taxid, seen_as in sorted(tax.unresolved.items()):
        print(f"    {taxid:>10s}  rank={tax.rank.get(taxid, '?'):<14s} "
              f"{seen_as or tax.name(taxid) or '?'}")


def report_tally(calls):
    """Per-profiler verdict counts, in VERDICTS order."""
    tally = {}
    for row in calls:
        counts = tally.setdefault(row["profiler"], dict.fromkeys(VERDICTS, 0))
        counts[row["verdict"]] += 1
    for profiler, counts in tally.items():
        print(f"  {profiler:12s} " + "  ".join(f"{v}={counts[v]}" for v in VERDICTS))


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True,
                    help="aMeta results dir (contains AUTHENTICATION/, KRAKENUNIQ/)")
    ap.add_argument("--truth", required=True,
                    help="per-read simulation truth: either a directory holding "
                         "<sample>/truth.tsv, or a path pattern containing "
                         "{sample}. This, not the config composition tables, is "
                         "the ground truth -- see the module docstring.")
    ap.add_argument("--taxonomy", default="",
                    help="directory holding NCBI nodes.dmp (and ideally names.dmp, "
                         "merged.dmp), or the path to nodes.dmp itself. Required: "
                         "taxa are joined by species taxid and there is no "
                         "name-based fallback. Default: auto-detect under "
                         "<results>/NGSLCA_DB, MALT_DB or KRAKENUNIQ_DB.")
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
    return ap.parse_args(argv)


def discover_targets(auth_dir, profilers_arg):
    """(profilers, samples): the sample universe is the union across profilers."""
    if profilers_arg:
        profilers = [p for p in profilers_arg.split(",") if p]
    else:
        profilers = sorted(d for d in os.listdir(auth_dir)
                           if os.path.isdir(os.path.join(auth_dir, d)))
    samples = set()
    for p in profilers:
        pdir = os.path.join(auth_dir, p)
        samples |= {s for s in os.listdir(pdir)
                    if os.path.isdir(os.path.join(pdir, s))}

    def natural(s):
        m = re.search(r"(\d+)$", s)
        return (int(m.group(1)) if m else 1 << 30, s)
    return profilers, sorted(samples, key=natural)


def summarise_profiler(profiler, args, auth_dir, samples, truth, host_keys,
                       detected, tax):
    """Render one profiler's heatmap and return its (metric rows, call rows)."""
    scores = load_profiler_scores(auth_dir, profiler, samples, tax)
    rows, label, smax = build_matrix(scores, truth, samples,
                                     args.include_host, host_keys, tax)
    render_heatmap(
        profiler, rows, label, samples, scores, truth, smax, detected,
        os.path.join(args.outdir, f"auth_summary_heatmap.{profiler}.pdf"),
        os.path.join(args.outdir, f"auth_summary_heatmap.{profiler}.png"))

    # keep integral targets integral: argparse applies type=float only to values
    # it actually parses, so otherwise the diff columns would switch between
    # "56" and "56.0" depending on whether --modern-score was passed
    perfect = {"ancient": args.ancient_score, "modern": args.modern_score}
    perfect = {k: int(v) if float(v).is_integer() else v
               for k, v in perfect.items()}
    per_sample = {
        s: sample_metrics(s, scores, truth, detected, perfect,
                          args.include_host, host_keys, auth_dir, profiler)
        for s in samples
    }
    return (metric_rows(profiler, samples, per_sample),
            call_rows(profiler, samples, scores, truth, detected, label,
                      args.include_host, host_keys))


def main(argv=None):
    args = parse_args(argv)

    auth_dir = os.path.join(args.results, "AUTHENTICATION")
    if not os.path.isdir(auth_dir):
        sys.exit(f"No AUTHENTICATION dir under {args.results}")

    profilers, samples = discover_targets(auth_dir, args.profilers)
    tax = open_taxonomy(args.taxonomy, args.results)
    truth, host_keys = load_truth(args.truth, samples, tax)
    detected = load_kraken_detected(args.results, samples, tax)

    os.makedirs(args.outdir, exist_ok=True)
    metrics, calls = [], []
    for profiler in profilers:
        m, c = summarise_profiler(profiler, args, auth_dir, samples, truth,
                                  host_keys, detected, tax)
        metrics += m
        calls += c

    report_unresolved(tax)

    out_calls = os.path.join(args.outdir, "auth_calls.tsv")
    write_tsv(out_calls, CALL_COLS, calls)
    print(f"\nPer-call detail -> {out_calls}")
    report_tally(calls)

    places = dict(METRIC_SPEC)
    out_metrics = os.path.join(args.outdir, "auth_metrics.tsv")
    write_tsv(out_metrics, METRIC_COLS, metrics, places)
    print(f"\nMetrics -> {out_metrics}")
    print_tsv(METRIC_COLS, metrics, places)


if __name__ == "__main__":
    main()
