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
  4. `auth_accession_conflicts.tsv`: reference sequences the two profilers
     assign to different taxa, with the false-positive reads each put on them.
     malt-build labels a reference by parsing taxon names out of its FASTA
     description (`-tn/--parseTaxonNames`, default true) as well as looking the
     accession up, and takes the LCA of every name it matches; ngsLCA only looks
     the accession up. A description that names a second organism -- a homolog
     it is "similar to", a cloning host, a lab in a clone name -- therefore
     pulls MALT to a different node, often a coarser rank, so the same reference
     can drive a species-level false positive for one profiler and a
     genus-level (unscored, invisible) call for the other. That makes the two
     profilers' false-positive taxon counts non-comparable at face value, and
     this table is the correction term.
  5. `auth_calls.tsv`: one row per (profiler, sample, taxon) with read support,
     score, breadth of coverage and a TP/FP/FN verdict -- the per-cell backing
     for both of the above. It doubles as a cache: breadth is the one column that
     costs a full pass over a `samtools depth -a` file per node, so on a rerun an
     existing auth_calls.tsv in --outdir is reused, per profiler, whenever every
     other column of it still matches what this run computes (see
     `cache_for_profiler`). Pass --no-cache to always remeasure.

Every result that leans on KrakenUniq's screen is *KrakenUniq-adjusted* by
default: a taxon the screen never surfaced was never offered to the profiler, so
it counts neither as a prediction nor as a miss and the numbers describe the
authentication step alone. `--no-kraken-adjust` switches the three summaries that
report a single adjusted figure -- the Jaccard curve, the criteria tables and the
ANI scatter's y axis -- to the unadjusted view, which charges the profiler for
the screening step in front of it as well, and writes them under `.unadjusted`
names so the two views can sit side by side in one outdir. The rest already
carries both readings and so does not move: auth_metrics.tsv has `diff_full`
beside `diff_kraken_adjusted`, auth_calls.tsv has the `krakenuniq_surfaced`
column and splits misses into FN_profiler / FN_krakenuniq, auth_ani.tsv has both
deficit columns, and the heatmap gives a never-surfaced cell its own style. That
also keeps auth_calls.tsv, and with it the breadth-of-coverage cache, identical
between the two modes.

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
import collections
import csv
import gzip
import importlib.util
import shutil
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
             "verdict", "breadth", "evenness"]

# The auth_calls.tsv columns measured off a node's per-base depth file, and how
# each is rendered. Both come from one pass, so both are cached together.
COVERAGE_COLS = (("breadth", "{:.6f}"), ("evenness", "{:.4f}"))

# A cached row is identified by these, and reusable only if every column that is
# not `breadth` still agrees -- i.e. the samples, their composition and what the
# profiler did with them are unchanged.
CACHE_KEY_COLS = ("profiler", "sample", "taxon")
COMPOSITION_COLS = [c for c in CALL_COLS
                    if c not in [name for name, _ in COVERAGE_COLS]]

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


# meta_bench/genomes, i.e. next to this script's parent -- the FASTAs the
# simulation drew from, which the ANI comparison measures the profilers'
# authentication references against.
DEFAULT_GENOMES = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "genomes")


# --------------------------------------------------------------------------- #
# Accession-level taxonomy conflicts between profilers
# --------------------------------------------------------------------------- #
def taxon_label(taxid, tax):
    """'Nonomuraea [genus] (83681)'.

    Rank is part of the label on purpose: the two sides routinely disagree
    about *rank* rather than about lineage, and a bare name hides that.
    """
    if not taxid:
        return ""
    return f"{tax.name(taxid) or '?'} [{tax.rank.get(taxid, 'no rank')}] ({taxid})"


def auth_bam_accessions(bam):
    """({accession: distinct reads}, {accession: taxid}) for one sorted.bam.

    Authentication BAMs name their references two ways, and the difference is
    the whole point of this table: MALT writes `<accession>|tax|<taxid>`, which
    is *its own DB's* label for that sequence, while the aligner+ngsLCA path
    writes a bare accession and leaves labelling to ngsLCA's -acc2tax. So one
    profiler's assignment arrives with the reads and the other has to be looked
    up (load_acc2taxid).
    """
    reads, label = {}, {}
    if not os.path.isfile(bam):
        return {}, label
    proc = subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, text=True)
    with proc:
        for line in proc.stdout:
            f = line.split("\t", 3)
            if len(f) < 3 or f[2] == "*":
                continue
            acc, _, tagged = f[2].partition("|tax|")
            if tagged:
                label[acc] = tagged.split("|")[0]
            reads.setdefault(acc, set()).add(f[0])
    return {a: len(r) for a, r in reads.items()}, label


def accession_fp_reads(auth_dir, profiler, calls, scores):
    """({accession: false-positive reads}, {accession: this profiler's taxid}).

    Only nodes this profiler called FP are counted: a node's read_list/sorted.bam
    under AUTHENTICATION is exactly the read support behind its score, so this
    breaks the read support of the *wrong* calls down by the reference sequence
    that attracted them.
    """
    reads, label = collections.Counter(), {}
    for row in calls:
        if row["profiler"] != profiler or row["verdict"] != "FP":
            continue
        rec = scores.get((row["sample"], row["_key"]))
        if not rec:
            continue
        per_acc, tags = auth_bam_accessions(
            os.path.join(auth_dir, profiler, row["sample"], rec["taxid"],
                         "sorted.bam"))
        for acc, n in per_acc.items():
            reads[acc] += n
        label.update(tags)
    return reads, label


def profiler_sam_path(results_dir, profiler, sample):
    """The canonical SAM for a profiler+sample, mirroring aMeta's auth_alignment_sam."""
    sub = "MALT" if profiler == "malt" else "ALIGNMENT"
    return os.path.join(results_dir, sub, f"{sample}.trimmed.sam.gz")


def sam_header_labels(path):
    """({accession: taxid}, n_refs) from a SAM's @SQ block.

    MALT stamps its own label into every reference name (`<acc>|tax|<taxid>`),
    so its whole DB labelling is readable straight off the header -- which sits
    entirely at the top of the file, ahead of the first alignment, so this reads
    a few MB rather than the whole 200 MB SAM. A reference MALT could not label
    at all carries no tag and is simply absent from the result. The aligner path
    writes bare accessions and yields nothing here, which is how the caller
    tells the two apart.
    """
    labels, refs = {}, 0
    if not os.path.isfile(path):
        return labels, refs
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith("@"):
                break
            if not line.startswith("@SQ"):
                continue
            refs += 1
            for field in line.rstrip("\n").split("\t"):
                if not field.startswith("SN:"):
                    continue
                acc, _, tagged = field[3:].partition("|tax|")
                if tagged:
                    labels[acc] = tagged.split("|")[0]
    return labels, refs


def load_acc2taxid(path, needed):
    """{accession.version: taxid} restricted to `needed`, streamed.

    Handles both layouts this pipeline keeps around: NCBI's
    nucl_gb.accession2taxid (header row, accession.version in column 2) and the
    two-column seqid2taxid map the reference FASTA is subset with. The first is
    ~10 GB, hence the streaming and the early exit.
    """
    found = {}
    if not path or not os.path.isfile(path) or not needed:
        return found
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        wide = len(head) >= 3
        if not head or head[0] != "accession":
            fh.seek(0)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if wide:
                if len(f) < 3:
                    continue
                acc, taxid = f[1], f[2]
            elif len(f) >= 2:
                acc, taxid = f[0], f[1]
            else:
                continue
            if acc in needed:
                found[acc] = taxid
                if len(found) == len(needed):
                    break
    return found


def accession_conflict_cols(profilers):
    a, b = profilers
    return ["accession", f"{a}_assignment", f"{b}_assignment",
            f"{a}_fp_reads", f"{b}_fp_reads"]


def accession_conflict_rows(profilers, labels, reads, tax):
    """One row per accession the two profilers label with different taxids."""
    a, b = profilers
    rows, unlabelled = [], 0
    for acc in sorted(set(reads[a]) | set(reads[b])):
        ta, tb = labels[a].get(acc), labels[b].get(acc)
        if not ta or not tb:
            unlabelled += 1
            continue
        if ta == tb:
            continue
        rows.append({
            "accession": acc,
            f"{a}_assignment": taxon_label(ta, tax),
            f"{b}_assignment": taxon_label(tb, tax),
            f"{a}_fp_reads": reads[a].get(acc, 0),
            f"{b}_fp_reads": reads[b].get(acc, 0),
        })
    rows.sort(key=lambda r: (-(r[f"{a}_fp_reads"] + r[f"{b}_fp_reads"]),
                             r["accession"]))
    return rows, unlabelled


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


def node_table(path):
    """{column: value} for a one-data-row MaltExtract table, or {} if unreadable.

    Handles both layouts: MaltExtract writes a rowname header ("Taxon" plus the
    columns), while the ngsLCA extractor omits it and its data row is one field
    longer than the header -- so the two cannot simply be zipped by a DictReader.
    """
    if not os.path.isfile(path):
        return {}
    with open(path) as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if any(r)]
    if len(rows) < 2:
        return {}
    header, data = rows[0], rows[1]
    vals = data[1:] if len(data) == len(header) + 1 else data
    return dict(zip(header, vals))


def read_node_reads(folder):
    """Aligned reads on the node's reference (TotalAlignmentsOnReference), or 0."""
    rd = os.path.join(folder, "MaltExtract_output", "default", "readDist")
    if not os.path.isdir(rd):
        return 0
    for fn in sorted(os.listdir(rd)):
        if fn.endswith("_alignmentDist.txt"):
            cell = node_table(os.path.join(rd, fn)).get(
                "TotalAlignmentsOnReference", 0)
            try:
                return int(float(cell))
            except (ValueError, TypeError):
                return 0
    return 0


# breadth and score.R's evenness statistic in one pass: tile the reference into
# 100 percentile bins, count the bases sitting in bins that are under 1% covered,
# and report that in bins' worth (score.R lines 29-41, `sum(V4<0.01)/step`).
AWK_COVERAGE = (
    'BEGIN{step=(hi-lo)/100}'
    '{t++; if($3>0) c++;'
    ' if(step>0){i=int($2/step); if(i<100){n[i]++; if($3>0) k[i]++}}}'
    'END{if(t==0) exit;'
    ' for(i=0;i<100;i++) if(n[i] && k[i]/n[i]<0.01) bad+=n[i];'
    ' printf "%.6f\\t%s", c/t, (step>0 ? sprintf("%.4f", bad/step) : "NA")}'
)


def last_line(path, tail_bytes=4096):
    """The last non-empty line, read from the end rather than by scanning."""
    with open(path, "rb") as fh:
        fh.seek(0, os.SEEK_END)
        fh.seek(max(0, fh.tell() - tail_bytes))
        lines = fh.read().decode(errors="replace").splitlines()
    for line in reversed(lines):
        if line.strip():
            return line
    return ""


def node_coverage(folder):
    """(breadth, evenness) for one node, from a single pass over its depth file.

    breadth is the fraction of reference positions with depth > 0 -- the "X% of
    genome covered" figure aMeta prints on its evenness plot. evenness is the
    statistic score.R thresholds at < 3 for its two evenness points. Both come off
    the same `samtools depth -a` file, which runs to hundreds of MB, so they are
    measured together; the tile width needs the first and last position, and those
    are read directly off the ends instead of by scanning. Either may be None.
    """
    path = os.path.join(folder, "breadth_of_coverage")
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None, None
    try:
        with open(path) as fh:
            first = fh.readline().split("\t")
        last = last_line(path).split("\t")
        lo, hi = int(first[1]), int(last[1])
    except (OSError, IndexError, ValueError):
        return None, None
    try:
        r = subprocess.run(["awk", "-v", f"lo={lo}", "-v", f"hi={hi}",
                            AWK_COVERAGE, path], capture_output=True, text=True)
    except OSError:
        return None, None
    fields = r.stdout.strip().split("\t")
    if len(fields) != 2:
        return None, None

    def parse(text):
        try:
            return None if text == "NA" else float(text)
        except ValueError:
            return None
    return parse(fields[0]), parse(fields[1])


def fmt_coverage(coverage, sample, key, index):
    """Render one coverage cell: blank if never measured, else NA or the value."""
    if (sample, key) not in coverage:
        return ""
    value = coverage[(sample, key)][index]
    return "NA" if value is None else COVERAGE_COLS[index][1].format(value)


def load_call_cache(path):
    """A previous auth_calls.tsv as profiler -> {(key cols): row}, or {} if unusable."""
    if not os.path.isfile(path):
        return {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c, _ in COVERAGE_COLS
                   if c not in (reader.fieldnames or [])]
        if missing:
            print(f"Coverage cache: ignoring {path} (no {', '.join(missing)} "
                  f"column -- written by an older version of this script)")
            return {}
        rows = list(reader)
    cache = {}
    for row in rows:
        key = tuple(row.get(c, "") for c in CACHE_KEY_COLS)
        cache.setdefault(row.get("profiler", ""), {})[key] = row
    print(f"Coverage cache: {path} ({len(rows)} row(s), {len(cache)} profiler(s))")
    return cache


def cache_for_profiler(cache, profiler, calls):
    """Cached coverage for `profiler` as (sample, taxon) -> (breadth, evenness).

    Coverage is measured off one `samtools depth` file per authenticated node, so
    a cached value only carries over when this run is looking at the same nodes
    of the same samples. Rather than guess which columns could move independently
    of the depth files, this demands that *every* other column match: the same
    (sample, taxon) rows, the same simulated read counts and fractions, the same
    scores, aligned reads and verdicts. Any disagreement discards the whole
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
              f"({len(complaints)} row(s)), remeasuring coverage:")
        for line in complaints[:5]:
            print(f"      {line}")
        if len(complaints) > 5:
            print(f"      ... and {len(complaints) - 5} more")
        return {}

    out = {}
    for (_p, sample, taxon), row in cached.items():
        texts = [(row.get(c) or "").strip() for c, _ in COVERAGE_COLS]
        if not any(texts):                 # never measured -- not a cache hit
            continue

        def parse(text):
            try:
                return None if text in ("", "NA") else float(text)
            except ValueError:
                return None
        out[(sample, taxon)] = tuple(parse(t) for t in texts)
    return out


def node_coverages(auth_dir, profiler, needed, scores, label, cached):
    """(sample, key) -> (breadth, evenness) for `needed`, reusing `cached`."""
    out, reused = {}, 0
    for sample, key in needed:
        hit = (sample, label.get(key, key))
        if hit in cached:
            out[(sample, key)] = cached[hit]
            reused += 1
            continue
        out[(sample, key)] = node_coverage(
            os.path.join(auth_dir, profiler, sample, scores[(sample, key)]["taxid"]))
    if reused:
        print(f"  [{profiler}] coverage: {reused} node(s) from cache, "
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
                   include_host, host_keys, coverage):
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
        cov = coverage.get((sample, key))
        if cov and cov[0] is not None:
            breadth[key in present].append(cov[0])

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

    Rows carry the join key under "_key" so callers can resolve breadth and a
    node's AUTHENTICATION directory; it is not in CALL_COLS, so it never
    reaches the TSV.
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
# Jaccard between the true ancient set and what a score threshold predicts
# --------------------------------------------------------------------------- #
def ancient_truth_set(sample, truth, host_keys, include_host):
    """The genomes simulated into `sample` as ancient."""
    return {k for k, info in truth.items()
            if info["category"] == "ancient" and sample in info["abund"]
            and (include_host or k not in host_keys)}


def jaccard_at_threshold(sample, scores, truth, host_keys, include_host, cut,
                         surfaced):
    """Jaccard of {scored >= cut} against the sample's true ancient set.

    The prediction is the raw score cut, with no appeal to ground truth -- unlike
    auth_metrics' jaccard_ancient, which intersects the found set with the ancient
    taxa before comparing and so cannot see a false positive that is a real but
    modern organism. Here every taxon clearing the cut counts as a claim that it
    is an ancient genome, which is what a score threshold actually asserts.

    `surfaced` is what KrakenUniq offered in this sample, and both sides are
    restricted to it -- the same adjustment diff_kraken_adjusted makes. It cuts
    in both directions: a node the profiler only holds because an earlier run's
    screen left its directory behind stops counting as a prediction, and a genome
    KrakenUniq never offered stops counting as a miss. What is left measures the
    authentication step rather than the screening step in front of it.

    Pass surfaced=None (--no-kraken-adjust) for the unadjusted variant, which
    compares against the sample's whole true ancient set and so scores the
    screening step and the authentication step together, as a user of the
    pipeline experiences them.
    """
    predicted = {k for (s, k), rec in scores.items()
                 if s == sample and rec["score"] >= cut
                 and (surfaced is None or k in surfaced)
                 and (include_host or k not in host_keys)}
    actual = ancient_truth_set(sample, truth, host_keys, include_host)
    return _jaccard(predicted, actual if surfaced is None else actual & surfaced)


def jaccard_threshold_rows(profilers, scores_by_profiler, samples, truth,
                           host_keys, include_host, cuts, detected, adjust=True):
    """One row per threshold, each profiler's mean-over-samples Jaccard."""
    rows = []
    for cut in cuts:
        row = {"threshold": cut}
        for p in profilers:
            row[p] = _mean([
                jaccard_at_threshold(s, scores_by_profiler[p], truth, host_keys,
                                     include_host, cut,
                                     detected.get(s, set()) if adjust else None)
                for s in samples])
        rows.append(row)
    return rows


def render_jaccard_curve(rows, profilers, out_pdf, out_png, n_samples,
                         adjusted=True):
    """Mean Jaccard against score threshold, one line per profiler."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    cmap = plt.get_cmap("tab10")
    markers = ["o", "s", "^", "D", "v"]

    for i, p in enumerate(profilers):
        xs = [r["threshold"] for r in rows]
        ys = [r[p] for r in rows]
        ax.plot(xs, ys, marker=markers[i % len(markers)], markersize=5,
                linewidth=1.6, color=cmap(i), label=p, zorder=3)

    ax.set_xticks([r["threshold"] for r in rows])
    ax.set_xlabel("authentication score threshold t  (predicted ancient = score ≥ t)")
    ax.set_ylabel(f"mean Jaccard over {n_samples} "
                  f"{'sample' if n_samples == 1 else 'samples'}")
    ax.set_title("Recovery of the true ancient genome set vs score threshold\n"
                 + ("(KrakenUniq-adjusted)" if adjusted else ""),
                 fontsize=11)
    ax.set_ylim(0, max(1.0, max((r[p] for r in rows for p in profilers
                                 if r[p] == r[p]), default=0) * 1.15))
    ax.grid(True, linewidth=0.4, color="0.88", zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.legend(fontsize=9, frameon=False, loc="best")

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Jaccard-vs-threshold curve -> {out_pdf}")


# --------------------------------------------------------------------------- #
# score.R's criteria, per node
# --------------------------------------------------------------------------- #
# The ten points score.R awards, and what ground truth says each should do. It is
# an *ancient* DNA authentication score, so the expectation depends on what the
# simulation put in the sample:
#
#   absent  -- nothing should fire: the organism is not there at all.
#   ancient -- everything should fire: it is present, damaged and short.
#   modern  -- the four damage-dependent criteria should NOT fire, and the rest
#              should. Failing them is the score working correctly, so counting
#              those failures as misses would invert the meaning of the table.
#              Read length is one of the four: it is a proxy for ancient
#              fragmentation, and gargammel's modern pass draws longer fragments
#              (empirical sizefreq, 83% under 100bp) than its ancient one
#              (lognormal, median ~42bp), so a modern genome cannot reach it.
#
# (name, points, expected of an ancient genome, expected of a modern one)
CRITERIA = (
    ("damage C>T_1",            1, True, False),
    ("damage G>A_20",           1, True, False),
    ("evenness of coverage",    2, True, True),
    ("edit distance (all)",     1, True, True),
    ("edit distance (ancient)", 1, True, False),
    ("PMD",                     1, True, False),
    ("read length",             1, True, False),
    ("depth of coverage",       1, True, True),
    ("ANI",                     1, True, True),
)
CRITERIA_EXPECTED = {"ancient": [c[2] for c in CRITERIA],
                     "modern": [c[3] for c in CRITERIA],
                     "absent": [False] * len(CRITERIA)}


def _num(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return float("nan")


def _decreasing(values):
    return all(values[i] > values[i + 1] for i in range(len(values) - 1))


def node_criteria(folder, rma6, evenness):
    """Which of the CRITERIA this node satisfies, as a list of bools.

    A transcription of score.R against the same files it reads. `evenness` is
    passed in because it is measured off the per-base depth file alongside
    breadth (see node_coverage) rather than re-read here.
    """
    m = os.path.join(folder, "MaltExtract_output")
    dam = node_table(f"{m}/default/damageMismatch/{rma6}_damageMismatch.txt")
    ed_all = node_table(f"{m}/default/editDistance/{rma6}_editDistance.txt")
    ed_anc = node_table(f"{m}/ancient/editDistance/{rma6}_editDistance.txt")
    dist = node_table(f"{m}/default/readDist/{rma6}_alignmentDist.txt")
    pid = node_table(f"{m}/default/percentIdentity/{rma6}_percentIdentity.txt")

    pmd_path = os.path.join(folder, "PMDscores.txt")
    pmd = False
    if os.path.isfile(pmd_path) and os.path.getsize(pmd_path):
        scores = [_num(l.split("\t")[3]) for l in open(pmd_path)
                  if len(l.split("\t")) > 3]
        pmd = bool(scores) and sum(v > 3 for v in scores) / len(scores) > 0.1

    lengths_path = os.path.join(folder, "read_length.txt")
    short = False
    if os.path.isfile(lengths_path):
        lengths = [_num(l) for l in open(lengths_path) if l.strip()]
        short = bool(lengths) and sum(v < 100 for v in lengths) / len(lengths) > 0.9

    a, b, c = (_num(pid.get(k)) for k in ("90", "95", "100"))
    total = a + b + c
    ani = total == total and total != 0 and \
        (a * 90 + b * 95 + c * 100) / total > 97

    return [
        round(_num(dam.get("C>T_1")), 4) > 0.05,
        round(_num(dam.get("G>A_20")), 4) > 0.05,
        evenness is not None and evenness < 3,
        _decreasing([_num(ed_all.get(str(i))) for i in range(5)]),
        _decreasing([_num(ed_anc.get(str(i))) for i in range(1, 5)]),
        pmd,
        short,
        _num(dist.get("TotalAlignmentsOnReference")) > 200,
        ani,
    ]


def criteria_counts(profiler, auth_dir, samples, scores, truth, coverage,
                    host_keys, include_host, detected, keep, adjust=True):
    """Per-criterion (FP, FN, expected-pass, expected-fail, skipped) for a profiler.

    Every node the profiler scored is a trial. Its ground-truth category fixes
    what each criterion should have done, and the node's own tables say what it
    did: a criterion that fired when it should not have is a false positive, one
    that stayed silent when it should have fired is a false negative.

    The trials are KrakenUniq-adjusted at both ends unless `adjust` is off
    (--no-kraken-adjust), which judges every node the profiler scored -- including
    any left behind by an earlier run's screen, see below. An in-sample taxon KrakenUniq
    never surfaced has no node, so it is never charged a false negative -- it was
    never the authentication step's to find. In the other direction a node whose
    taxon is absent from *this* run's taxID.species is skipped: aMeta creates node
    directories with `mkdir -p` and never prunes them, and its checkpoints glob the
    directory listing rather than the species file, so a results directory reused
    across runs keeps authenticating (and scoring) whatever an earlier KrakenUniq
    screen turned up. Judging those would charge the profiler for a leftover
    directory.

    A node scoring 0 failed everything -- score.R only ever adds -- so it is
    counted without touching the filesystem, which is also why coverage is only
    ever measured for the nodes that scored.
    """
    n = len(CRITERIA)
    fp, fn, exp_pass, exp_fail = [0] * n, [0] * n, [0] * n, [0] * n
    skipped = 0
    for (sample, key), rec in scores.items():
        if not include_host and key in host_keys:
            continue
        if adjust and key not in detected.get(sample, set()):
            skipped += 1
            continue
        category = (truth[key]["category"] if is_present(truth, key, sample)
                    else "absent")
        if not keep(key, category):
            continue
        expected = CRITERIA_EXPECTED[category]
        if rec["score"] > 0:
            cov = coverage.get((sample, key))
            got = node_criteria(
                os.path.join(auth_dir, profiler, sample, rec["taxid"]),
                f"{sample}.trimmed.rma6", cov[1] if cov else None)
        else:
            got = [False] * n
        for i, (want, have) in enumerate(zip(expected, got)):
            exp_pass[i] += want
            exp_fail[i] += not want
            fp[i] += have and not want
            fn[i] += want and not have
    return fp, fn, exp_pass, exp_fail, skipped


def criteria_rows(profilers, auth_dir, samples, scores_by_profiler, truth,
                  coverage_by_profiler, host_keys, include_host, detected, keep,
                  adjust=True):
    """(rows, skipped): one row per criterion, each profiler's FP/FN side by side."""
    per_profiler = {
        p: criteria_counts(p, auth_dir, samples, scores_by_profiler[p], truth,
                           coverage_by_profiler[p], host_keys, include_host,
                           detected, keep, adjust)
        for p in profilers
    }
    rows = []
    for i, (name, points, _anc, _mod) in enumerate(CRITERIA):
        row = {"criterion": name, "points": points}
        for p in profilers:
            fp, fn, exp_pass, exp_fail, _skipped = per_profiler[p]
            row[f"{p}_FP"] = fp[i]
            row[f"{p}_FN"] = fn[i]
            row[f"{p}_expected_pass"] = exp_pass[i]
            row[f"{p}_expected_fail"] = exp_fail[i]
        rows.append(row)
    return rows, {p: per_profiler[p][4] for p in profilers}


# FP reuses the heatmap's "not in sample" red, so a false positive reads the same
# across both figures; FN gets its own amber rather than the heatmap's blue, which
# already means "modern".
FP_COLOR = CAT_COLOR["absent"]
FN_COLOR = (0.90, 0.55, 0.10)


# The criteria tables rendered on every run: (filename suffix, title suffix, the
# ground-truth categories judged).
#
# The "ancient" view drops modern genomes but keeps the absent ones. Restricting
# to ancient *nodes* alone would produce a dead table: ancient expects all nine
# criteria to fire, so "fired when it should not have" cannot happen and the whole
# FP column would be structurally 0. Absent nodes are what make false positives
# observable, so they stay; dropping only the modern genomes is what makes this
# the ancient-DNA view, and it is the same reasoning that keeps modern taxa off
# the ANI scatter.
CRITERIA_VIEWS = (
    ("", "", ("ancient", "modern", "absent")),
    ("ancient", " (ancient genomes)", ("ancient", "absent")),
)


def criteria_filter(categories, taxids):
    """A keep(key, category) predicate for one criteria view."""
    def keep(key, category):
        return category in categories and (not taxids or key in taxids)
    return keep


def expected_label(index):
    """Which ground-truth categories should fire criterion `index`."""
    _name, _pts, ancient, modern = CRITERIA[index]
    return "ancient + modern" if ancient and modern else "ancient only"


def render_criteria_table(rows, profilers, out_pdf, out_png, judged, skipped,
                          title_suffix=""):
    """Draw the FP/FN table as a figure, cells shaded by count."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    widths = [2.5, 0.42, 1.35] + [0.62, 0.62] * len(profilers)
    edges = [0.0]
    for w in widths:
        edges.append(edges[-1] + w)
    span = edges[-1]

    nbody = len(rows)
    nrow = nbody + 3                       # 2 header rows + body + total
    rh = 0.30                              # row height, same units as widths
    fig_w, fig_h = span + 0.3, nrow * rh + 1.15
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, span)
    ax.set_ylim(0, nrow * rh)
    ax.axis("off")

    def y(row):                            # row 0 is the top header
        return (nrow - 1 - row) * rh

    def cell(col, row, text, *, face=None, weight="normal", align="center",
             color="black", size=8, ncols=1):
        x0, x1 = edges[col], edges[col + ncols]
        if face is not None:
            ax.add_patch(Rectangle((x0, y(row)), x1 - x0, rh, facecolor=face,
                                   edgecolor="white", lw=0.8, zorder=1))
        pad = 0.06
        tx = {"center": (x0 + x1) / 2, "left": x0 + pad, "right": x1 - pad}[align]
        ax.text(tx, y(row) + rh / 2, text, ha=align, va="center", fontsize=size,
                color=color, fontweight=weight, zorder=2)

    header = (0.90, 0.90, 0.92)
    cell(0, 0, "", face=header)
    cell(1, 0, "", face=header)
    cell(2, 0, "should fire for", face=header, weight="bold", size=7.5)
    for i, p in enumerate(profilers):
        cell(3 + 2 * i, 0, p, face=header, weight="bold", ncols=2)
    for col, text in ((0, "criterion"), (1, "pts")):
        cell(col, 1, text, face=header, weight="bold",
             align="left" if col == 0 else "center")
    cell(2, 1, "", face=header)
    for i, _p in enumerate(profilers):
        cell(3 + 2 * i, 1, "FP", face=header, weight="bold", color=FP_COLOR)
        cell(4 + 2 * i, 1, "FN", face=header, weight="bold", color=FN_COLOR)

    max_fp = max((r[f"{p}_FP"] for r in rows for p in profilers), default=1) or 1
    max_fn = max((r[f"{p}_FN"] for r in rows for p in profilers), default=1) or 1
    for j, row in enumerate(rows):
        r = j + 2
        stripe = (1, 1, 1) if j % 2 else (0.965, 0.965, 0.975)
        cell(0, r, row["criterion"], face=stripe, align="left")
        cell(1, r, str(row["points"]), face=stripe)
        cell(2, r, expected_label(j), face=stripe, size=7, color="0.35")
        for i, p in enumerate(profilers):
            for col, value, base, top in ((3 + 2 * i, row[f"{p}_FP"], FP_COLOR, max_fp),
                                          (4 + 2 * i, row[f"{p}_FN"], FN_COLOR, max_fn)):
                face = blend(base, value / top) if value else stripe
                cell(col, r, str(value), face=face,
                     color="white" if luminance(face) < 0.55 else "black")

    r = nbody + 2
    cell(0, r, "total", face=header, weight="bold", align="left")
    cell(1, r, str(sum(x["points"] for x in rows)), face=header, weight="bold")
    cell(2, r, "", face=header)
    for i, p in enumerate(profilers):
        cell(3 + 2 * i, r, str(sum(x[f"{p}_FP"] for x in rows)), face=header,
             weight="bold")
        cell(4 + 2 * i, r, str(sum(x[f"{p}_FN"] for x in rows)), face=header,
             weight="bold")

    ax.set_title("Authentication criteria results" + title_suffix,
                 fontsize=11, pad=12)
    note = ("FP = criterion fired where ground truth says it should not; "
            "FN = vice versa")
            # + "   ".join(f"{p}: {judged[p]} node(s) judged"
            #              + (f", {skipped[p]} skipped (not in this run's "
            #                 f"KrakenUniq output)" if skipped[p] else "")
            #              for p in profilers))
    fig.text(0.5, 0.012, note, ha="center", va="bottom", fontsize=7, color="0.35")

    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Criteria table -> {out_pdf}")


def print_criteria_table(rows, profilers, judged, skipped):
    """The FP/FN table, profilers across the top and criteria down the side."""
    width = max(len(r["criterion"]) for r in rows) + 2
    cell = 6
    head1 = " " * (width + 5) + "|" + "|".join(
        p.center(cell * 2 + 1) for p in profilers) + "|"
    head2 = f"{'criterion':<{width}}{'pts':>4} |" + "|".join(
        f"{'FP':>{cell}}{'FN':>{cell + 1}}" for _ in profilers) + "|"
    print("\n" + head1)
    print(head2)
    print("-" * len(head2))
    for r in rows:
        line = f"{r['criterion']:<{width}}{r['points']:>4} |"
        for p in profilers:
            line += f"{r[f'{p}_FP']:>{cell}}{r[f'{p}_FN']:>{cell + 1}}|"
        print(line)
    print("-" * len(head2))
    totals = f"{'total':<{width}}{sum(r['points'] for r in rows):>4} |"
    for p in profilers:
        totals += (f"{sum(r[f'{p}_FP'] for r in rows):>{cell}}"
                   f"{sum(r[f'{p}_FN'] for r in rows):>{cell + 1}}|")
    print(totals)
    for p in profilers:
        extra = (f", {skipped[p]} skipped (not in this run's KrakenUniq output)"
                 if skipped.get(p) else "")
        print(f"  {p}: {judged[p]} node(s) judged{extra}")


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

# The scatter's y axis: a row key `ani_rows` fills with whichever of the two
# deficits this run plots, so the label placer and the connectors below never
# have to know which mode they are in. It is absent from ANI_COLS, so it stays
# out of the TSV -- which carries both deficits either way.
#
# The KrakenUniq-adjusted deficit is the default because it is the one that says
# something about authentication: a taxon KrakenUniq never surfaced was never
# offered to the profiler, so the whole authentication chain never ran for it and
# its full deficit measures the screening step instead of the scoring. Under
# --no-kraken-adjust the full deficit is plotted instead and the x axis keeps its
# meaning, but a point can then sit high because KrakenUniq never offered the
# taxon rather than because its reference is divergent.
Y_METRIC = "_plotted_diff"


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
             host_keys, genomes_dir, perfect, include_host, tax, min_af,
             adjust=True):
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

    `adjust` picks which deficit the scatter plots (see Y_METRIC); both are in
    every row regardless. It does not change the taxon universe, which stays the
    taxa KrakenUniq surfaced at least once: a taxon it never surfaced anywhere
    has no node directory and so no authentication reference to compare, which
    would exclude it in either mode.
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
                Y_METRIC: adjusted if adjust else full,
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


def render_ani_scatter(rows, profilers, out_pdf, out_png, n_samples,
                       adjusted=True):
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
    ax.set_ylabel(("mean KrakenUniq-adjusted " if adjusted else "mean ")
                  + "|perfect − score| over "
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
    ap.add_argument("--criteria-taxids", default="",
                    help="comma-separated species taxids to restrict the "
                         "criteria tables to (default: every judged node). The "
                         "taxid is the species-level join key, i.e. the row "
                         "label on the ANI scatter and the `taxid` column of "
                         "auth_ani.tsv.")
    ap.add_argument("--no-ani", action="store_true",
                    help="skip the skani ANI comparison and its scatter plot")
    ap.add_argument("--min-aligned-frac", type=float, default=0.15,
                    help="minimum fraction of both genomes skani must align for "
                         "an ANI to reach the scatter (default 0.15); below it "
                         "the top reference is a fragment rather than a genome "
                         "and its ANI says nothing. Such rows stay in "
                         "auth_ani.tsv marked plotted=no.")
    ap.add_argument("--acc2taxid", default="",
                    help="accession -> taxid mapping used to label references "
                         "for the profiler that does not tag its BAM RNAMEs "
                         "(the aligner+ngsLCA path). Pass the file ngsLCA was "
                         "run with -- the config's malt_accession2taxid -- for "
                         "its exact labels, at the cost of streaming a ~10 GB "
                         "file (minutes). Default: "
                         "<results>/REFERENCE_DB/seqid2taxid.project.map, which "
                         "is the same mapping for almost every accession, is "
                         "read in seconds, but is not the file ngsLCA reads.")
    ap.add_argument("--no-accession-conflicts", action="store_true",
                    help="skip auth_accession_conflicts.tsv (it costs one "
                         "samtools pass per false-positive node)")
    ap.add_argument("--no-cache", action="store_true",
                    help="always remeasure breadth of coverage instead of "
                         "reusing an auth_calls.tsv already in --outdir")
    ap.add_argument("--no-kraken-adjust", action="store_true",
                    help="Report the kraken-unajusted plots. Default (adjusted) is "
                         "the authentication-only view -- a taxon not included in the screen "
                         " counts neither as a prediction nor as a miss. ")
    ap.add_argument("--include-host", action="store_true",
                    help="keep host taxa (e.g. Homo sapiens); default drops them")
    return ap.parse_args(argv)


def adjust_tag(args):
    """Filename tag that keeps the two KrakenUniq views out of each other's way."""
    return ".unadjusted" if args.no_kraken_adjust else ""


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
    # first, coverage-less; only measured nodes (scored, and plotted at all) are
    # worth a depth pass -- the host's depth file is the largest of the lot. A
    # node that scored 0 failed every criterion by construction, so it needs no
    # measurement for either the metrics or the criteria table.
    calls = call_rows(profiler, samples, scores, truth, detected, label,
                      args.include_host, host_keys)
    needed = {(s, k) for (s, k), rec in scores.items()
              if rec["score"] > 0 and (args.include_host or k not in host_keys)}
    coverage = node_coverages(auth_dir, profiler, needed, scores, label,
                              cache_for_profiler(cache, profiler, calls))
    for row in calls:
        # kept on the row (render_row reads only the explicit columns, so it
        # reaches neither the TSV nor the cache comparison) because later
        # consumers need the node's taxid to find its AUTHENTICATION directory
        key = row["_key"]
        for index, (name, _fmt) in enumerate(COVERAGE_COLS):
            row[name] = fmt_coverage(coverage, row["sample"], key, index)

    perfect = perfect_scores(args)
    per_sample = {
        s: sample_metrics(s, scores, truth, detected, perfect,
                          args.include_host, host_keys, coverage)
        for s in samples
    }
    return metric_rows(profiler, samples, per_sample), calls, scores, coverage


def summarise_jaccard(args, profilers, scores_by_profiler, samples, truth,
                      host_keys, detected):
    """Write and draw the Jaccard-vs-score-threshold curve."""
    # Thresholds start at 1 and run to the perfect ancient score. ">= 0" would
    # predict every node the profiler screened, including the ones it explicitly
    # scored 0 and rejected, which is not a threshold anyone would set; ">= 10"
    # is the other end, "only a perfect score counts".
    adjust = not args.no_kraken_adjust
    tag = adjust_tag(args)
    cuts = list(range(1, int(args.ancient_score) + 1))
    rows = jaccard_threshold_rows(profilers, scores_by_profiler, samples, truth,
                                  host_keys, args.include_host, cuts, detected,
                                  adjust)
    out_tsv = os.path.join(args.outdir, f"auth_jaccard_threshold{tag}.tsv")
    write_tsv(out_tsv, ["threshold"] + list(profilers), rows,
              {p: 4 for p in profilers})
    print(f"\nJaccard vs score threshold ({'KrakenUniq-adjusted' if adjust else 'unadjusted'})"
          f" -> {out_tsv}")
    render_jaccard_curve(
        rows, profilers,
        os.path.join(args.outdir, f"auth_jaccard_threshold{tag}.pdf"),
        os.path.join(args.outdir, f"auth_jaccard_threshold{tag}.png"),
        len(samples), adjust)
    print_tsv(["threshold"] + list(profilers), rows, {p: 4 for p in profilers})


def summarise_criteria(args, profilers, auth_dir, samples, scores_by_profiler,
                       truth, coverage_by_profiler, host_keys, detected):
    """Write and draw every CRITERIA_VIEWS table, honouring --criteria-taxids."""
    taxids = {t.strip() for t in args.criteria_taxids.split(",") if t.strip()}
    if taxids:
        print(f"\nCriteria tables restricted to {len(taxids)} taxid(s): "
              + ", ".join(sorted(taxids)))
    cols = ["criterion", "points"] + [
        f"{p}_{suffix}" for p in profilers
        for suffix in ("FP", "FN", "expected_pass", "expected_fail")]

    for slug, title_suffix, categories in CRITERIA_VIEWS:
        rows, skipped = criteria_rows(
            profilers, auth_dir, samples, scores_by_profiler, truth,
            coverage_by_profiler, host_keys, args.include_host, detected,
            criteria_filter(categories, taxids), not args.no_kraken_adjust)
        judged = {p: rows[0][f"{p}_expected_pass"] + rows[0][f"{p}_expected_fail"]
                  for p in profilers}
        if not any(judged.values()):
            print(f"\nCriteria table{title_suffix or ' (all genomes)'}: "
                  "no node matched the filter, skipping")
            continue

        # <name>.<variant>.<ext>, the same shape as auth_summary_heatmap.malt.pdf
        tag = (f".{slug}" if slug else "") + adjust_tag(args)
        out_tsv = os.path.join(args.outdir, f"auth_criteria{tag}.tsv")
        write_tsv(out_tsv, cols, rows)
        print(f"\nPer-criterion FP/FN{title_suffix} -> {out_tsv}")
        render_criteria_table(
            rows, profilers,
            os.path.join(args.outdir, f"auth_criteria_table{tag}.pdf"),
            os.path.join(args.outdir, f"auth_criteria_table{tag}.png"),
            judged, skipped, title_suffix)
        print_criteria_table(rows, profilers, judged, skipped)


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

    adjust = not args.no_kraken_adjust
    rows, skipped = ani_rows(
        profilers, scores_by_profiler, auth_dir, samples, truth, detected,
        host_keys, args.genomes, perfect_scores(args), args.include_host, tax,
        args.min_aligned_frac, adjust)

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

    # the table is the same in both modes -- it carries both deficits -- so only
    # the scatter, which plots one of them, takes the .unadjusted name
    tag = adjust_tag(args)
    out_tsv = os.path.join(args.outdir, "auth_ani.tsv")
    write_tsv(out_tsv, ANI_COLS, rows, ANI_PLACES)
    print(f"  {len(rows)} row(s) -> {out_tsv}")
    render_ani_scatter(rows, profilers,
                       os.path.join(args.outdir, f"auth_ani_scatter{tag}.pdf"),
                       os.path.join(args.outdir, f"auth_ani_scatter{tag}.png"),
                       len(samples), adjust)


def profiler_reference_labels(results_dir, profiler, samples):
    """{accession: taxid} as *this profiler* labels its references, or {}.

    Read from the @SQ block of the profiler's own SAM, not from the per-node
    authentication BAMs: those hold only the one reference a node was
    authenticated against, which would restrict the comparison to references
    both profilers happened to pick as a top hit and so hide exactly the
    disagreements that make one profiler call a species where the other calls a
    genus. Returns {} for a profiler that writes bare accessions -- its
    labelling lives in the acc2taxid file instead.
    """
    labels = {}
    for i, sample in enumerate(samples):
        found, refs = sam_header_labels(
            profiler_sam_path(results_dir, profiler, sample))
        labels.update(found)
        # an untagged profiler stays untagged for every sample; stop after the
        # first rather than re-reading a million-line header nine more times
        if i == 0 and refs and not found:
            break
    return labels


def summarise_accession_conflicts(args, profilers, auth_dir, calls,
                                  scores_by_profiler, tax):
    """Write auth_accession_conflicts.tsv: references the profilers label apart.

    malt-build resolves each reference by *parsing taxon names out of the FASTA
    description* (`-tn/--parseTaxonNames`, default true) alongside its -a2t
    accession lookup, and takes the LCA of everything it matches; ngsLCA does a
    plain accession -> taxid lookup and nothing else. Where a description names
    a second organism -- a homolog it is "similar to", a cloning host, or a lab
    in a clone name ("Cancer Center" -> the crab genus Cancer) -- the two land
    on different nodes, usually at different ranks. The same reference then
    feeds a species-level false positive on one side and a genus-level, and so
    unscored and invisible, call on the other, which is why the profilers'
    false-positive taxon counts are not comparable at face value.

    Columns: the accession, what each profiler calls it, and the false-positive
    reads each profiler actually put on it (0 meaning that profiler aligned
    false-positive reads elsewhere, not that the reference went unused).
    """
    if len(profilers) != 2:
        print(f"\nAccession conflicts skipped: needs exactly 2 profilers, "
              f"got {len(profilers)} ({', '.join(profilers) or 'none'})")
        return
    if shutil.which("samtools") is None:
        print("\nAccession conflicts skipped: samtools is not on PATH")
        return

    print("\nAccession-level taxonomy conflicts:")
    samples = sorted({row["sample"] for row in calls})
    reads, labels = {}, {}
    for profiler in profilers:
        reads[profiler], _tags = accession_fp_reads(
            auth_dir, profiler, calls, scores_by_profiler[profiler])
        labels[profiler] = profiler_reference_labels(
            args.results, profiler, samples)
        src = (f"{len(labels[profiler])} from its SAM @SQ tags"
               if labels[profiler] else "untagged, needs acc2taxid")
        print(f"  {profiler}: {len(reads[profiler])} reference(s) carrying "
              f"false-positive reads; labels {src}")

    universe = set(reads[profilers[0]]) | set(reads[profilers[1]])
    untagged = [p for p in profilers if not labels[p]]
    if untagged:
        path = args.acc2taxid or os.path.join(
            args.results, "REFERENCE_DB", "seqid2taxid.project.map")
        found = load_acc2taxid(path, universe)
        print(f"  acc2taxid: {len(found)}/{len(universe)} accession(s) resolved "
              f"from {path}")
        if not args.acc2taxid:
            print("    (that is the reference-subset map, not the file ngsLCA "
                  "reads -- pass --acc2taxid <config's malt_accession2taxid> "
                  "for its exact labels)")
        for profiler in untagged:
            labels[profiler] = found

    rows, unlabelled = accession_conflict_rows(profilers, labels, reads, tax)
    cols = accession_conflict_cols(profilers)
    out = os.path.join(args.outdir, "auth_accession_conflicts.tsv")
    write_tsv(out, cols, rows)
    print(f"  {len(rows)} conflicting reference(s) of {len(universe)} carrying "
          f"false-positive reads -> {out}")
    if unlabelled:
        print(f"  {unlabelled} reference(s) not compared: only one side had a "
              "label (MALT leaves a reference untagged when it can resolve no "
              "name for it)")
    if rows:
        a, b = profilers
        print("  top conflicts by false-positive read support:")
        for row in rows[:8]:
            print(f"    {row['accession']:<16} "
                  f"{a}={row[f'{a}_assignment']} [{row[f'{a}_fp_reads']} reads]"
                  f"  |  {b}={row[f'{b}_assignment']} "
                  f"[{row[f'{b}_fp_reads']} reads]")


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

    metrics, calls = [], []
    scores_by_profiler, coverage_by_profiler = {}, {}
    for profiler in profilers:
        m, c, s, cov = summarise_profiler(profiler, args, auth_dir, samples,
                                          truth, host_keys, detected, tax, cache)
        metrics += m
        calls += c
        scores_by_profiler[profiler] = s
        coverage_by_profiler[profiler] = cov

    summarise_criteria(args, profilers, auth_dir, samples, scores_by_profiler,
                       truth, coverage_by_profiler, host_keys, detected)
    summarise_jaccard(args, profilers, scores_by_profiler, samples, truth,
                      host_keys, detected)

    if not args.no_ani:
        summarise_ani(args, profilers, scores_by_profiler, auth_dir, samples,
                      truth, detected, host_keys, tax)

    if not args.no_accession_conflicts:
        summarise_accession_conflicts(args, profilers, auth_dir, calls,
                                      scores_by_profiler, tax)

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
