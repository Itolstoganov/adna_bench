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
  3. `auth_ani.tsv` + `auth_ani_scatter.{pdf,png}`: skani ANI between the
     reference a taxon was authenticated against and the genome the simulation
     actually drew its reads from, against that taxon's mean KrakenUniq-adjusted
     score deficit (the TSV carries the unadjusted one too). One point per
     (profiler, taxon) over the *ancient* taxa KrakenUniq surfaced at least once
     -- a modern genome's deficit is set by the damage criteria it cannot reach,
     not by its reference, so it is tabulated but not plotted.
     This is the diagnostic for score deficits that are a property of the
     reference database rather than of the profiler: several score.R criteria
     (edit distance, ANI, evenness) compare reads to the *top reference*, which
     is only the simulated genome when that genome is in the aMeta reference DB.
  4. `auth_calls.tsv`: one row per (profiler, sample, taxon) with read support,
     score, breadth of coverage and a TP/FP/FN verdict -- the per-cell backing
     for both of the above. It doubles as a cache: breadth is the one column that
     costs a full pass over a `samtools depth -a` file per node, so on a rerun an
     existing auth_calls.tsv in --outdir is reused, per profiler, whenever every
     other column of it still matches what this run computes (see
     `cache_for_profiler`). Pass --no-cache to always remeasure.

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
import importlib.util
import math
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
             "verdict", "breadth"]

# A cached row is identified by these, and reusable only if every column that is
# not `breadth` still agrees -- i.e. the samples, their composition and what the
# profiler did with them are unchanged.
CACHE_KEY_COLS = ("profiler", "sample", "taxon")
COMPOSITION_COLS = [c for c in CALL_COLS if c != "breadth"]

# meta_bench/genomes, i.e. next to this script's parent -- the FASTAs the
# simulation drew from, which the ANI comparison measures the profilers'
# authentication references against.
DEFAULT_GENOMES = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "genomes")


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
      key -> {category, display, genome, abund:{sample: frac}, reads:{sample: n}}
    and `frac` is the taxon's share of that sample's reads. `genome` keeps the
    raw truth.tsv label, which is the basename of the FASTA gargammel drew from
    (`Sorangium_cellulosum` -> genomes/Sorangium_cellulosum.fa) and so is what
    locates the simulated genome again for the ANI comparison; `display` is the
    same thing prettified for the tables and plots.
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
                "genome": row.taxon,
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


def fmt_breadth(breadths, sample, key):
    """Render a breadth cell: blank when the node was never measured, else NA/value."""
    if (sample, key) not in breadths:
        return ""
    value = breadths[(sample, key)]
    return "NA" if value is None else f"{value:.6f}"


def load_call_cache(path):
    """A previous auth_calls.tsv as profiler -> {(key cols): row}, or {} if unusable."""
    if not os.path.isfile(path):
        return {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "breadth" not in reader.fieldnames:
            print(f"Breadth cache: ignoring {path} (no breadth column -- written "
                  f"by an older version of this script)")
            return {}
        rows = list(reader)
    cache = {}
    for row in rows:
        key = tuple(row.get(c, "") for c in CACHE_KEY_COLS)
        cache.setdefault(row.get("profiler", ""), {})[key] = row
    print(f"Breadth cache: {path} ({len(rows)} row(s), {len(cache)} profiler(s))")
    return cache


def cache_for_profiler(cache, profiler, calls):
    """Cached breadths for `profiler` as (sample, taxon) -> float|None, or {}.

    Breadth is measured off one `samtools depth` file per authenticated node, so
    a cached value only carries over when this run is looking at the same nodes
    of the same samples. Rather than guess which columns could move independently
    of the depth files, this demands that *every* non-breadth column match: the
    same (sample, taxon) rows, the same simulated read counts and fractions, the
    same scores, aligned reads and verdicts. Any disagreement discards the whole
    profiler's cache -- a partial reuse would be reporting numbers from two
    different runs in one table.
    """
    cached = cache.get(profiler)
    if not cached:
        return {}

    complaints, fresh_keys = [], set()
    for row in calls:
        key = tuple(str(row[c]) for c in CACHE_KEY_COLS)
        fresh_keys.add(key)
        old = cached.get(key)
        if old is None:
            complaints.append(f"{key[1]}/{key[2]}: absent from the cached table")
        elif (render_row(row, COMPOSITION_COLS, {})
              != [old.get(c, "") for c in COMPOSITION_COLS]):
            complaints.append(f"{key[1]}/{key[2]}: composition or score changed")
    complaints += [f"{k[1]}/{k[2]}: cached but no longer present"
                   for k in sorted(set(cached) - fresh_keys)]

    if complaints:
        print(f"  [{profiler}] cached calls do not match this run "
              f"({len(complaints)} row(s)), remeasuring breadth:")
        for line in complaints[:5]:
            print(f"      {line}")
        if len(complaints) > 5:
            print(f"      ... and {len(complaints) - 5} more")
        return {}

    out = {}
    for (_p, sample, taxon), row in cached.items():
        text = (row.get("breadth") or "").strip()
        if not text:                       # never measured -- not a cache hit
            continue
        try:
            out[(sample, taxon)] = None if text == "NA" else float(text)
        except ValueError:
            continue
    return out


def node_breadths(auth_dir, profiler, needed, scores, label, cached):
    """(sample, key) -> breadth for the nodes in `needed`, reusing `cached`."""
    out, reused = {}, 0
    for sample, key in needed:
        hit = (sample, label.get(key, key))
        if hit in cached:
            out[(sample, key)] = cached[hit]
            reused += 1
            continue
        out[(sample, key)] = node_breadth(
            os.path.join(auth_dir, profiler, sample, scores[(sample, key)]["taxid"]))
    if reused:
        print(f"  [{profiler}] breadth: {reused} node(s) from cache, "
              f"{len(needed) - reused} measured")
    return out


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
                   include_host, host_keys, breadths):
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
        b = breadths.get((sample, key))
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
    """One row per (sample, taxon) that is either simulated or authenticated.

    Rows carry the join key under "_key" for the caller to resolve breadth with;
    it is not in CALL_COLS, so it never reaches the TSV.
    """
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
                "_key": key,
            })
    return rows


# --------------------------------------------------------------------------- #
# ANI: authentication reference vs simulated genome
# --------------------------------------------------------------------------- #
ANI_COLS = ["profiler", "taxid", "taxon", "category", "auth_ref", "sim_genome",
            "ani", "aligned_frac_sim", "aligned_frac_auth",
            "mean_diff_kraken_adjusted", "mean_diff_full", "mean_truth_frac",
            "abund_intensity", "samples_present", "samples_scored", "plotted"]
ANI_PLACES = {"ani": 2, "aligned_frac_sim": 3, "aligned_frac_auth": 3,
              "mean_diff_kraken_adjusted": 2, "mean_diff_full": 2,
              "mean_truth_frac": 5, "abund_intensity": 3}

# The scatter's y axis. The KrakenUniq-adjusted deficit is the one that says
# something about authentication: a taxon KrakenUniq never surfaced was never
# offered to the profiler, so the whole authentication chain never ran for it and
# its full deficit measures the screening step instead of the scoring.
Y_METRIC = "mean_diff_kraken_adjusted"


def read_contigs(path):
    """Yield each FASTA record's sequence as one bytes object, for pyskani."""
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield b"".join(seq)
                name, seq = line[1:].strip(), []
            elif name is not None:
                seq.append(line.strip().encode())
    if name is not None:
        yield b"".join(seq)


def skani_ani(auth_fasta, sim_fasta, cache):
    """(ANI%, frac of sim aligned, frac of auth aligned), or None if no hit.

    Uses the pyskani bindings rather than the skani CLI so the measurement needs
    no external process and no temporary sketch files. skani screens out pairs it
    cannot align (below roughly 80% ANI, or too little shared sequence) by
    returning no hit at all -- that is reported as None rather than as a number,
    since "too divergent to align" is not the same as a low ANI estimate.
    """
    hit = cache.get((auth_fasta, sim_fasta), False)
    if hit is not False:
        return hit
    import pyskani
    db = pyskani.Database()
    db.sketch("auth", *read_contigs(auth_fasta))
    hits = db.query("sim", *read_contigs(sim_fasta))
    # a self-consistent pair yields one hit; keep the best-aligned one otherwise
    best = max(hits, key=lambda h: h.query_fraction, default=None)
    out = None if best is None else (best.identity * 100.0,
                                     best.query_fraction, best.reference_fraction)
    cache[(auth_fasta, sim_fasta)] = out
    return out


def abundance_intensity(key, truth, samples, smax):
    """The heatmap's colour intensity for one taxon, averaged over its samples.

    `build_matrix` normalises each heatmap cell against the most abundant taxon
    in *that* sample, so a cell's intensity reads as "how much of this sample is
    this organism, relative to its biggest member". Averaging that same ratio
    over the samples the taxon was simulated into carries the reading to a
    scatter point, where there is one mark per taxon rather than one per cell.
    """
    abund = truth[key]["abund"]
    return _mean([abund[s] / smax[s] for s in samples if s in abund]) or 0.0


def sample_abundance_max(truth, samples, host_keys, include_host):
    """Per-sample maximum ground-truth abundance: the intensity normaliser.

    Matches the `smax` build_matrix computes for the heatmaps -- the host dwarfs
    every microbe, so including it when it is not drawn would flatten the range.
    """
    smax = {}
    for sample in samples:
        vals = [info["abund"].get(sample, 0.0) for k, info in truth.items()
                if include_host or k not in host_keys]
        smax[sample] = max(vals, default=0.0) or 1.0
    return smax


def auth_reference_fasta(auth_dir, profiler, samples, key, scores):
    """(ref_id, path) of the reference this taxon was authenticated against.

    Breadth_Of_Coverage writes the top reference of each node next to it as
    <ref_id>.fasta, so the exact sequence score.R judged the reads against is
    already on disk. Nodes with no reads are dropped before that rule runs, so
    the first sample that has the file wins.
    """
    for sample in samples:
        rec = scores.get((sample, key))
        if rec is None:
            continue
        node = os.path.join(auth_dir, profiler, sample, rec["taxid"])
        name_list = os.path.join(node, "name_list.txt")
        if not os.path.isfile(name_list):
            continue
        with open(name_list) as fh:
            ref_id = fh.readline().strip()
        path = os.path.join(node, f"{ref_id}.fasta")
        if ref_id and os.path.isfile(path):
            return ref_id, path
    return None, None


def taxon_mean_diff(key, samples, scores, truth, perfect, detected):
    """(KrakenUniq-adjusted, full) mean |perfect - score| for one taxon.

    Both run over every sample: samples where the taxon is absent contribute
    their own deficit (0 for a correct rejection, the score itself for a false
    positive), which is why the mean is not taken over the in-sample cells only.
    That makes each value one taxon's contribution to the matching
    auth_metrics.tsv column, scaled by the sample count -- summing over *all*
    taxa would reproduce that metric, though auth_ani.tsv holds only the ones
    with both a reference and a simulated genome and so will not add up to it.

    The adjusted variant drops the cells where the taxon was simulated but
    KrakenUniq never surfaced it -- the same rule diff_kraken_adjusted uses.
    Those are misses of the screening step, which happens before any reference is
    chosen, so counting them against reference divergence would be reading a
    KrakenUniq result off an axis that knows nothing about KrakenUniq.
    """
    full = adjusted = 0
    for sample in samples:
        in_sample = is_present(truth, key, sample)
        want = perfect[truth[key]["category"]] if in_sample else 0
        rec = scores.get((sample, key))
        d = abs(want - (rec["score"] if rec else 0))
        full += d
        if not (in_sample and key not in detected.get(sample, set())):
            adjusted += d
    n = len(samples)
    if not n:
        return float("nan"), float("nan")
    return adjusted / n, full / n


def ani_rows(profilers, scores_by_profiler, auth_dir, samples, truth, detected,
             host_keys, genomes_dir, perfect, include_host, tax, min_af):
    """One row per (profiler, taxon) for every taxon KrakenUniq ever surfaced.

    A taxon only earns a row when both ends of the comparison exist: the
    authentication reference (so the profiler actually scored it somewhere) and
    the simulated genome (so there is something to compare it against). A
    KrakenUniq hit on a taxon the simulation never contained -- a false positive
    -- has no simulated genome by definition and is reported as skipped rather
    than plotted against an invented x.

    `plotted` marks the rows the scatter uses; the rest stay in the TSV with
    their numbers intact. A row is held back when

      * the taxon is modern in the sample. score.R is an ancient-DNA score, so a
        modern genome forfeits both damage criteria, the PMD criterion and (the
        modern gargammel pass drawing longer fragments) the read-length one --
        four of ten points, none of which the reference can influence. Its
        deficit measures the damage model, and reading it off a divergence axis
        would be reading an answer to a different question.
      * skani aligned less than `min_af` of either genome, because the top
        reference is not always a whole genome: MALT picks e.g. the 1.1 kb
        KJ477415.1 for Lactococcus lactis, over which skani reports a perfectly
        true and perfectly meaningless 100% ANI.
      * skani found no hit at all.
    """
    universe = set()
    for keys in detected.values():
        universe |= keys
    if not include_host:
        universe -= host_keys

    rows, skipped, cache = [], {}, {}
    smax = sample_abundance_max(truth, samples, host_keys, include_host)
    for key in sorted(universe, key=lambda k: (truth.get(k, {}).get("display")
                                               or tax.name(k) or k).lower()):
        info = truth.get(key)
        if info is None:
            skipped.setdefault("not simulated (false positive)", []).append(key)
            continue
        sim_fasta = os.path.join(genomes_dir, f"{info['genome']}.fa")
        if not os.path.isfile(sim_fasta):
            skipped.setdefault(f"no genome under {genomes_dir}", []).append(key)
            continue
        for profiler in profilers:
            scores = scores_by_profiler[profiler]
            ref_id, auth_fasta = auth_reference_fasta(
                auth_dir, profiler, samples, key, scores)
            if auth_fasta is None:
                skipped.setdefault("never authenticated (no top reference)",
                                   []).append(key)
                continue
            ani = skani_ani(auth_fasta, sim_fasta, cache)
            adjusted, full = taxon_mean_diff(key, samples, scores, truth,
                                             perfect, detected)
            if info["category"] != "ancient":
                reason = ("modern in the sample: four criteria are unreachable "
                          "without deamination, so its deficit is set by the "
                          "damage model rather than by the reference")
            elif ani is None:
                reason = "no skani hit (too divergent to align)"
            elif min(ani[1], ani[2]) < min_af:
                reason = (f"skani aligned < {min_af:.0%} of a genome (top "
                          "reference is not a whole genome)")
            else:
                reason = None
            if reason:
                skipped.setdefault(reason, []).append(key)
            rows.append({
                "profiler": profiler,
                "taxid": key,
                "taxon": info["display"],
                "category": info["category"],
                "auth_ref": ref_id,
                "sim_genome": os.path.basename(sim_fasta),
                "ani": float("nan") if ani is None else ani[0],
                "aligned_frac_sim": float("nan") if ani is None else ani[1],
                "aligned_frac_auth": float("nan") if ani is None else ani[2],
                "mean_diff_kraken_adjusted": adjusted,
                "mean_diff_full": full,
                "mean_truth_frac": _mean([info["abund"][s] for s in samples
                                          if s in info["abund"]]),
                "abund_intensity": abundance_intensity(key, truth, samples, smax),
                "samples_present": sum(is_present(truth, key, s) for s in samples),
                "samples_scored": sum((s, key) in scores for s in samples),
                "plotted": "no" if reason else "yes",
            })
    return rows, skipped


# An identical pair has no logarithm, so the divergence axis needs a floor. skani
# resolves differences well below this, so a point sitting on it really is "the
# same sequence" rather than a value that was rounded away.
DIVERGENCE_FLOOR = 1e-6
LABEL_FONTSIZE = 8


def divergence(row):
    """1 - ANI as a fraction: 99.9% ANI -> 1e-3. Floored so it can be logged."""
    return max(1.0 - row["ani"] / 100.0, DIVERGENCE_FLOOR)


def abbreviate_binomial(name):
    """'Sorangium cellulosum' -> 'S. cellulosum'; anything else is left alone."""
    parts = name.split()
    if len(parts) >= 2 and len(parts[0]) > 1:
        return f"{parts[0][0]}. " + " ".join(parts[1:])
    return name


def point_label(row):
    """'56 S. cellulosum': the taxid to join against the tables, the abbreviated
    binomial to read the plot without them."""
    return f"{row['taxid']} {abbreviate_binomial(row['taxon'])}".strip()


# A label further than this from its point (in points) gets a leader line: close
# enough and the pairing is obvious, far enough and it has to be drawn.
LEADER_THRESHOLD_PT = 12


def place_labels(ax, fig, groups, anchor_of, points, label_of):
    """Annotate each taxon once, greedily avoiding overlaps between labels.

    The points cluster hard, so a fixed offset would stack labels on top of each
    other. Each label tries offsets around its anchor in order of increasing
    displacement and takes the first whose (approximate) box misses everything
    placed so far; labels on the right of the axis try their left-hand offsets
    first so they stay inside. Anything that ends up displaced gets a leader line
    back to its point, which is what keeps the pairing readable once the placer
    has had to push a label past its neighbours.
    """
    fig.canvas.draw()                       # transData is only valid once drawn
    px_per_pt = fig.dpi / 72.0
    (x0_ax, y0_ax), (x1_ax, y1_ax) = ax.transAxes.transform([(0, 0), (1, 1)])

    # seed the occupied boxes with the markers themselves, so a label steps
    # around points as well as around other labels -- a name lying across a
    # marker is exactly the ambiguity the leader lines exist to prevent
    r = 5.5 * px_per_pt
    placed = []
    for x, y in points:
        mx, my = ax.transData.transform((x, y))
        placed.append((mx - r, my - r, mx + r, my + r))

    def overlaps(a, b):
        return not (a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1])

    # tallest first: the outliers are the labels most worth keeping legible
    for _taxid, group in sorted(groups.items(),
                                key=lambda kv: -max(r[Y_METRIC]
                                                    for r in kv[1])):
        anchor = anchor_of(group)
        label = label_of(anchor)
        px, py = ax.transData.transform((divergence(anchor),
                                         anchor[Y_METRIC]))
        w = 0.62 * LABEL_FONTSIZE * len(label) * px_per_pt
        h = 1.4 * LABEL_FONTSIZE * px_per_pt
        right_half = px > (x0_ax + x1_ax) / 2
        sides = ["right", "left"] if right_half else ["left", "right"]
        # nearest offsets first, so a label only travels as far as it must
        offsets = sorted(((dxm, dy) for dxm in (8, 22, 38)
                          for dy in (0, -12, 12, -24, 24, -36, 36)),
                         key=lambda o: math.hypot(*o))
        candidates = [(dxm if s == "left" else -dxm, dy, s)
                      for dxm, dy in offsets for s in sides]
        chosen = fallback = None
        for dx, dy, side in candidates:
            left = px + dx * px_per_pt - (w if side == "right" else 0)
            cy = py + dy * px_per_pt
            box = (left, cy - h / 2, left + w, cy + h / 2)
            # a label that runs off the axes is worse than one that overlaps
            if not (x0_ax <= box[0] and box[2] <= x1_ax
                    and y0_ax <= box[1] and box[3] <= y1_ax):
                continue
            if fallback is None:
                fallback = (dx, dy, side, box)
            if not any(overlaps(box, b) for b in placed):
                chosen = (dx, dy, side, box)
                break
        # every offset collided (or none fit): keep it in the axes and let it
        # overlap rather than dropping the label or letting it escape
        dx, dy, side, box = chosen or fallback or (8, 0, "left", None)
        if box is not None:
            placed.append(box)
        # a darkened marker colour, so the label also says which point is its own
        base = CAT_COLOR.get(anchor["category"], (0.25, 0.25, 0.25))
        ink = tuple(c * 0.62 for c in base)
        leader = None
        if math.hypot(dx, dy) > LEADER_THRESHOLD_PT:
            # leave from the label edge that faces the point, not its centre
            leader = dict(arrowstyle="-", linewidth=0.5, color=ink + (0.65,),
                          shrinkA=0, shrinkB=4,
                          relpos=(0.0 if side == "left" else 1.0, 0.5))
        ax.annotate(label, (divergence(anchor), anchor[Y_METRIC]),
                    textcoords="offset points", xytext=(dx, dy),
                    ha=side, va="center", fontsize=LABEL_FONTSIZE,
                    color=ink, zorder=4, arrowprops=leader)


def render_ani_scatter(rows, profilers, out_pdf, out_png, n_samples):
    """Score deficit (y) against reference divergence (x), one shape per profiler.

    ANI piles up between 98% and 100%, where a linear axis is unreadable, so x is
    the divergence 1 - ANI on a log scale -- each decade is a tenfold difference
    in how far the authentication reference sits from the genome the reads
    actually came from. A top axis carries the equivalent ANI. The two profilers'
    points for one taxon are joined by a thin line, so a taxon they disagree
    about (which is about the profiler) is distinguishable at a glance from one
    they agree on (which is about the reference).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plotted = [r for r in rows if r["plotted"] == "yes"]
    if not plotted:
        print("  no taxon has both an authentication reference and a simulated "
              "genome with a usable skani hit, skipping the ANI scatter")
        return

    markers = ["o", "s", "^", "D", "v"]
    shape = {p: markers[i % len(markers)] for i, p in enumerate(profilers)}
    fig, ax = plt.subplots(figsize=(9.5, 6.0))

    # keyed by taxid, not display name: the connector has to join the profilers'
    # views of one node, and two taxids can carry the same binomial. Labels are
    # "<taxid> <abbreviated binomial>" -- the taxid to join against auth_ani.tsv
    # and the heatmaps, the abbreviated name so the figure reads on its own
    # without a full binomial's width crowding the dense band.
    groups = {}
    for row in plotted:
        groups.setdefault(row["taxid"], []).append(row)

    for group in groups.values():
        if len(group) > 1:
            pair = sorted(group, key=divergence)
            ax.plot([divergence(r) for r in pair], [r[Y_METRIC] for r in pair],
                    color="0.72", linewidth=0.8, zorder=2)
    for row in plotted:
        base = CAT_COLOR.get(row["category"], (0.4, 0.4, 0.4))
        # same white -> category-colour ramp the heatmap cells use, so a pale
        # point here and a pale cell there mean the same thing. The edge is the
        # saturated colour rather than white: at low abundance the fill is nearly
        # background, and a white edge would leave nothing to see.
        ax.scatter(divergence(row), row[Y_METRIC],
                   marker=shape[row["profiler"]], s=58,
                   facecolor=blend(base, row["abund_intensity"]),
                   edgecolor=tuple(c * 0.62 for c in base),
                   linewidth=0.9, zorder=3)

    ax.set_xscale("log")
    xs = [divergence(r) for r in plotted]
    lo, hi = min(xs) / 2.5, max(xs) * 2.5
    ax.set_xlim(lo, hi)
    # Decade ticks carrying both readings: a secondary axis would inherit the log
    # locator and place its ANI ticks on the transformed values, which is nonsense.
    decades = [10.0 ** e for e in range(math.floor(math.log10(lo)),
                                        math.ceil(math.log10(hi)) + 1)]
    decades = [d for d in decades if lo <= d <= hi]
    ax.set_xticks(decades)
    ax.set_xticklabels([f"$10^{{{round(math.log10(d))}}}$\n{100 * (1 - d):g}% ANI"
                        for d in decades], fontsize=8)
    ax.set_xlabel("divergence 1 − ANI between the authentication reference and "
                  "the simulated genome (log scale)")
    ax.set_ylabel("mean KrakenUniq-adjusted |perfect − score| over "
                  f"{n_samples} {'sample' if n_samples == 1 else 'samples'}")
    ax.set_title("Authentication score deficit vs reference divergence "
                 "(ancient genomes)", fontsize=11)
    ax.grid(True, which="both", linewidth=0.4, color="0.88", zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    def swatch(face, edge, marker="o"):
        return Line2D([], [], linestyle="", marker=marker, markersize=7,
                      markerfacecolor=face, markeredgecolor=edge,
                      markeredgewidth=0.9, color="none")

    prof = [Line2D([], [], linestyle="", marker=shape[p], color="0.35",
                   markersize=7, label=p) for p in profilers]
    ax.add_artist(ax.legend(handles=prof, title="profiler", fontsize=8,
                            title_fontsize=8, frameon=False, loc="upper left"))

    # Anchor the ramp on real points rather than on 0/0.5/1: the taxa here occupy
    # a middle band of the scale (no organism is its sample's largest in every
    # sample, so averaging pulls everything inward), and a legend spanning the
    # full ramp would advertise a contrast the figure does not contain. Labelling
    # the swatches with the abundance they stand for makes it quantitative.
    ranked = sorted(plotted, key=lambda r: r["abund_intensity"])
    picks, seen = [], set()
    for row in (ranked[0], ranked[len(ranked) // 2], ranked[-1]):
        if row["taxid"] not in seen:
            seen.add(row["taxid"])
            picks.append(row)
    ramp, labels = [], []
    for row in picks:
        base = CAT_COLOR.get(row["category"], (0.4, 0.4, 0.4))
        ramp.append(swatch(blend(base, row["abund_intensity"]),
                           tuple(c * 0.62 for c in base)))
        labels.append(f"{row['mean_truth_frac'] * 100:.1f}% of a sample")
    ax.legend(handles=ramp, labels=labels, fontsize=8, title_fontsize=8,
              title="mean ground-truth abundance",
              frameon=False, loc="center left", handletextpad=0.6,
              labelspacing=0.35)

    # labels last: they need the final axis limits to know where the edges are
    place_labels(ax, fig, groups,
                 lambda g: max(g, key=lambda r: r[Y_METRIC]),
                 [(divergence(r), r[Y_METRIC]) for r in plotted],
                 point_label)

    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ANI scatter ({len(plotted)} point(s)) -> {out_pdf}")


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
    ap.add_argument("--modern-score", type=float, default=5,
                    help="perfect score for a modern genome (default 5, what one "
                         "can earn on merit: without deamination the two damage "
                         "criteria, the PMD criterion and -- because gargammel's "
                         "modern pass draws longer fragments -- the read-length "
                         "criterion are all unreachable. Use 6 for the empirical "
                         "ceiling, whose sixth point is the ancient edit-distance "
                         "criterion passing by chance on a handful of reads)")
    ap.add_argument("--genomes", default=DEFAULT_GENOMES,
                    help="directory of the FASTAs the simulation drew from, "
                         "named <taxon>.fa after the truth.tsv taxon column. "
                         "Used as the ANI scatter's reference point. "
                         f"Default: {DEFAULT_GENOMES}")
    ap.add_argument("--no-ani", action="store_true",
                    help="skip the skani ANI comparison and its scatter plot")
    ap.add_argument("--min-aligned-frac", type=float, default=0.15,
                    help="minimum fraction of both genomes skani must align for "
                         "an ANI to reach the scatter (default 0.15); below it "
                         "the top reference is a fragment rather than a genome "
                         "and its ANI says nothing. Such rows stay in "
                         "auth_ani.tsv marked plotted=no.")
    ap.add_argument("--no-cache", action="store_true",
                    help="always remeasure breadth of coverage instead of "
                         "reusing an auth_calls.tsv already in --outdir")
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


def perfect_scores(args):
    """The score a flawless authentication earns, per ground-truth category.

    Kept integral where possible: argparse applies type=float only to values it
    actually parses, so otherwise the diff columns would switch between "56" and
    "56.0" depending on whether --modern-score was passed.
    """
    perfect = {"ancient": args.ancient_score, "modern": args.modern_score}
    return {k: int(v) if float(v).is_integer() else v for k, v in perfect.items()}


def summarise_profiler(profiler, args, auth_dir, samples, truth, host_keys,
                       detected, tax, cache):
    """Render one profiler's heatmap, returning (metric rows, call rows, scores)."""
    scores = load_profiler_scores(auth_dir, profiler, samples, tax)
    rows, label, smax = build_matrix(scores, truth, samples,
                                     args.include_host, host_keys, tax)
    render_heatmap(
        profiler, rows, label, samples, scores, truth, smax, detected,
        os.path.join(args.outdir, f"auth_summary_heatmap.{profiler}.pdf"),
        os.path.join(args.outdir, f"auth_summary_heatmap.{profiler}.png"))

    # the calls are what a cached table is checked against, so they are built
    # first, breadth-less; only measured nodes (scored, and plotted at all) are
    # worth a depth pass -- the host's depth file is the largest of the lot.
    calls = call_rows(profiler, samples, scores, truth, detected, label,
                      args.include_host, host_keys)
    needed = {(s, k) for (s, k), rec in scores.items()
              if rec["score"] > 0 and (args.include_host or k not in host_keys)}
    breadths = node_breadths(auth_dir, profiler, needed, scores, label,
                             cache_for_profiler(cache, profiler, calls))
    for row in calls:
        row["breadth"] = fmt_breadth(breadths, row["sample"], row.pop("_key"))

    perfect = perfect_scores(args)
    per_sample = {
        s: sample_metrics(s, scores, truth, detected, perfect,
                          args.include_host, host_keys, breadths)
        for s in samples
    }
    return metric_rows(profiler, samples, per_sample), calls, scores


def summarise_ani(args, profilers, scores_by_profiler, auth_dir, samples, truth,
                  detected, host_keys, tax):
    """Write auth_ani.tsv and the ANI-vs-score-deficit scatter.

    Isolated from the rest of the run: pyskani is an extra dependency the other
    outputs do not need, so a missing binding costs the scatter and nothing else.
    """
    if importlib.util.find_spec("pyskani") is None:
        print("\nANI scatter skipped: pyskani is not installed in this "
              "environment.\n  pip install pyskani   (or rerun with --no-ani "
              "to stop asking)")
        return

    print(f"\nANI vs simulated genomes ({args.genomes}):")
    if not os.path.isdir(args.genomes):
        print(f"  no such directory, skipping -- pass --genomes")
        return

    rows, skipped = ani_rows(
        profilers, scores_by_profiler, auth_dir, samples, truth, detected,
        host_keys, args.genomes, perfect_scores(args), args.include_host, tax,
        args.min_aligned_frac)

    for reason, keys in sorted(skipped.items()):
        names = sorted({truth.get(k, {}).get("display") or tax.name(k) or k
                        for k in keys})
        print(f"  {len(names)} taxon(s) not plotted, {reason}: "
              + ", ".join(names[:6]) + (", ..." if len(names) > 6 else ""))

    # a taxon plotted for only some profilers gets no connector; say so, rather
    # than leaving a lone point looking like a rule about the x axis
    seen = {}
    for row in rows:
        if row["plotted"] == "yes":
            seen.setdefault(row["taxon"], []).append(row["profiler"])
    lone = {t: ps for t, ps in seen.items() if len(ps) < len(profilers)}
    if lone:
        print(f"  {len(lone)} taxon(s) plotted for only some profilers, so they "
              "carry no connector (the other row is listed above): "
              + ", ".join(f"{t} ({'+'.join(ps)} only)"
                          for t, ps in sorted(lone.items())))

    out_tsv = os.path.join(args.outdir, "auth_ani.tsv")
    write_tsv(out_tsv, ANI_COLS, rows, ANI_PLACES)
    print(f"  {len(rows)} row(s) -> {out_tsv}")
    render_ani_scatter(rows, profilers,
                       os.path.join(args.outdir, "auth_ani_scatter.pdf"),
                       os.path.join(args.outdir, "auth_ani_scatter.png"),
                       len(samples))


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
    out_calls = os.path.join(args.outdir, "auth_calls.tsv")
    cache = {} if args.no_cache else load_call_cache(out_calls)

    metrics, calls, scores_by_profiler = [], [], {}
    for profiler in profilers:
        m, c, s = summarise_profiler(profiler, args, auth_dir, samples, truth,
                                     host_keys, detected, tax, cache)
        metrics += m
        calls += c
        scores_by_profiler[profiler] = s

    if not args.no_ani:
        summarise_ani(args, profilers, scores_by_profiler, auth_dir, samples,
                      truth, detected, host_keys, tax)

    report_unresolved(tax)

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
