#!/usr/bin/env python3
"""Runtime and memory summary for an aMeta run, from Snakemake's benchmark files.

Every rule in this fork carries a `benchmark:` directive, so `benchmarks/<RULE>/
<sample>.benchmark.txt` already holds wall time, peak RSS and CPU time for each
job. This tool aggregates them into a per-stage comparison of the taxonomic
profilers, and adds the phase breakdown the tools print into their own logs --
which is the only way to see inside a rule.

Three tables:

  1. Per rule: wall time, peak RSS and CPU time over the samples that ran, plus
     the effective parallelism (cpu_time / wall) actually achieved. Threads are
     measured this way rather than read from the .smk, so the number reflects
     the run instead of the declaration.

  2. Per profiler: alignment and taxonomic assignment separately where the
     pipeline separates them, and their total. MALT does both in ONE process,
     so it has one rule and one benchmark file and cannot be split here -- see
     table 3 for the closest available breakdown.

  3. Phases, parsed out of the tools' own logs. This is where MALT's cost
     actually shows up: reloading the database index into a fresh JVM, once per
     sample, dwarfs the alignment itself.

Wall times are only comparable across rules when the jobs ran serially. Check
that before quoting a ratio: if Snakemake ran several samples at once, each
job's wall time includes the contention. Peak RSS for a JVM tool (MALT) is
bounded by its -Xmx heap and is an upper bound on what it needs, not a
requirement.

Aligner + ngsLCA profilers write to benchmarks/{ALIGNMENT,NGSLCA}/<aligner>/;
runs from before that layout (flat ALIGNMENT/, NGSLCA/) are read as strobealign.
Rules missing from the run (e.g. KrakenUniq reused from an earlier run) can be
taken from --extra-benchmarks; they are marked with '*'.

Usage:
    python bench_summary.py --rundir meta_bench/runs/sim_dataset
    python bench_summary.py --rundir <dir> --outdir <dir>   # also write TSVs + plots
    python bench_summary.py --rundir <dir> --extra-benchmarks <dir>/old_benchmarks
"""

import argparse
import csv
import glob
import os
import re
import statistics as stats
import sys

# (profiler, stage, rule) for the rules that make up a taxonomic profiler. A new
# backend needs one line here and nothing else.
PROFILER_STAGES = (
    ("malt", "alignment + LCA", "MALT"),
    ("strobealign", "alignment", "ALIGNMENT/strobealign"),
    ("strobealign", "taxonomic assignment", "NGSLCA/strobealign"),
    ("bowtie2", "alignment", "ALIGNMENT/bowtie2"),
    ("bowtie2", "taxonomic assignment", "NGSLCA/bowtie2"),
    ("bwa_aln", "alignment", "ALIGNMENT/bwa_aln"),
    ("bwa_aln", "taxonomic assignment", "NGSLCA/bwa_aln"),
    ("bwa_mem", "alignment", "ALIGNMENT/bwa_mem"),
    ("bwa_mem", "taxonomic assignment", "NGSLCA/bwa_mem"),
)
# Flat layout of single-aligner runs from before per-aligner directories.
LEGACY_RULES = {"ALIGNMENT/strobealign": "ALIGNMENT", "NGSLCA/strobealign": "NGSLCA"}
# Paid once per database/binary, not per sample.
SETUP_RULES = ("BUILD_MALT_DB", "COMPILE_STROBEALIGN", "REFERENCE_SUBSET",
               "STROBEALIGN_INDEX", "BOWTIE2_INDEX", "BWA_INDEX")
# Rules whose benchmark is not under benchmarks/<RULE>/, as globs relative to
# benchmarks/. Bowtie2_Index writes next to the FASTA it indexes; only the
# project reference's index is the profiler's (bowtie2_db is the pathogen screen).
BENCHMARK_GLOBS = {
    "BOWTIE2_INDEX": "../results/REFERENCE_DB/*_BOWTIE2_BUILD.benchmark.txt",
}
# The one-off index each profiler's alignment step reads (bwa aln and mem share one).
INDEX_RULES = {"malt": "BUILD_MALT_DB", "strobealign": "STROBEALIGN_INDEX",
               "bowtie2": "BOWTIE2_INDEX", "bwa_aln": "BWA_INDEX",
               "bwa_mem": "BWA_INDEX"}
# Upstream of every profiler, so it differentiates none of them -- but it sets
# the pipeline's real memory ceiling, so it is worth showing.
SHARED_RULES = ("KRAKENUNIQ", "FILTER_KRAKENUNIQ_OUTPUT", "CUTADAPT_ADAPTER_TRIMMING")

RULE_COLS = ["rule", "n", "wall_mean_s", "wall_min_s", "wall_max_s", "wall_total_s",
             "max_rss_mean_gb", "max_rss_max_gb", "cpu_time_mean_s", "eff_threads"]
PHASE_COLS = ["tool", "phase", "mean_s", "share_pct", "n", "note"]


def read_benchmarks(bench_dir, rule):
    """Rows of every benchmark file for `rule`. A repeat() benchmark has several."""
    rows = []
    pattern = BENCHMARK_GLOBS.get(rule, os.path.join(rule, "*.txt"))
    for path in sorted(glob.glob(os.path.join(bench_dir, pattern))):
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                row["_job"] = os.path.basename(path).replace(".benchmark.txt", "")
                rows.append(row)
    return rows


def _num(row, key, default=float("nan")):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return default


def summarise_rule(bench_dirs, rule):
    """Aggregate one rule's benchmark rows from the first of `bench_dirs` that
    has them, or None if it never ran. Rows from a fallback dir get a '*'."""
    rows = []
    for i, bench_dir in enumerate(bench_dirs):
        for name in (rule, LEGACY_RULES.get(rule)):
            if name and not rows:
                rows = read_benchmarks(bench_dir, name)
        if rows:
            break
    if not rows:
        return None
    wall = [_num(r, "s") for r in rows]
    rss = [_num(r, "max_rss") / 1024 for r in rows]          # Snakemake writes MB
    cpu = [_num(r, "cpu_time") for r in rows]
    total_wall = sum(wall)
    by_job = {}
    for r, w in zip(rows, wall):
        by_job.setdefault(r["_job"], []).append(w)
    return {
        "rule": rule + ("*" if i else ""),
        "n": len(rows),
        "wall_mean_s": stats.mean(wall),
        "wall_min_s": min(wall),
        "wall_max_s": max(wall),
        "wall_total_s": total_wall,
        "max_rss_mean_gb": stats.mean(rss),
        "max_rss_max_gb": max(rss),
        "cpu_time_mean_s": stats.mean(cpu),
        # what the run actually achieved, not what the rule asked for
        "eff_threads": (sum(cpu) / total_wall) if total_wall else float("nan"),
        # per sample, for plotting a profiler's stages summed sample by sample
        "wall_by_job": {job: stats.mean(ws) for job, ws in by_job.items()},
    }


# --------------------------------------------------------------------------- #
# Phase breakdown, parsed from the tools' own logs
# --------------------------------------------------------------------------- #
def _secs_after(lines, index, window=5):
    """The first '(12.3s)' timing within `window` lines after `index`."""
    for line in lines[index + 1:index + 1 + window]:
        m = re.search(r"\((\d+(?:\.\d+)?)s\)", line)
        if m:
            return float(m.group(1))
    return None


def malt_log_phases(path):
    """{phase: seconds} for one malt-run log.

    malt-run times the pieces it considers interesting and prints one 'Total
    time' at the end; LCA binning is NOT among the timed pieces, so it can only
    be reported lumped into the remainder together with JVM startup and the
    ncbi.map/ncbi.tre load.
    """
    lines = open(path, errors="replace").read().splitlines()
    out = {"db_load": 0.0, "align": 0.0, "write_sam": 0.0, "total": 0.0, "peak_gb": 0.0}
    for i, line in enumerate(lines):
        if line.startswith("Reading file:") and "maltDB.dat" in line:
            secs = _secs_after(lines, i, window=3)
            if secs is not None:
                out["db_load"] += secs
        elif line.startswith("+++++ Aligning file"):
            out["align"] = _secs_after(lines, i) or 0.0
        elif line.startswith("Copying from temporary file"):
            out["write_sam"] = _secs_after(lines, i, window=3) or 0.0
        elif line.startswith("Total time:"):
            m = re.search(r"([\d.]+)s", line)
            if m:
                out["total"] = float(m.group(1))
        elif line.startswith("Peak memory:"):
            m = re.search(r"([\d.]+)\s+of", line)
            if m:
                out["peak_gb"] = float(m.group(1))
    out["other"] = max(0.0, out["total"] - out["db_load"] - out["align"] - out["write_sam"])
    return out


def strobealign_log_phases(path):
    """{phase: seconds} for one strobealign log.

    Mind the units strobealign mixes: 'indexing' and 'mapping' are wall clock,
    while the finer breakdown underneath is summed across threads.
    """
    txt = open(path, errors="replace").read()

    def grab(label):
        m = re.search(re.escape(label) + r":\s*([\d.]+)\s*s", txt)
        return float(m.group(1)) if m else 0.0

    return {
        # indexing on the fly, or loading an index written by Strobealign_Index
        "index": grab("Total time indexing") or grab("Total time reading index"),
        "prebuilt": "Total time reading index" in txt,
        "map": grab("Total time mapping"),
        "hits_cpu": grab("Total time finding hits (non-rescue mode)"),
        "chain_cpu": grab("Total time chaining (non-rescue mode)"),
        "extend_cpu": grab("Total time extending and pairing seeds"),
    }


def ngslca_log_walltime(path):
    """ngsLCA's self-reported walltime, as a cross-check on the benchmark file."""
    m = re.search(r"walltime used\s*=\s*([\d.]+)\s*sec",
                  open(path, errors="replace").read())
    return float(m.group(1)) if m else None


def _logs(rundir, *patterns):
    """Logs matching the first of `patterns` that matches anything."""
    for pattern in patterns:
        paths = sorted(glob.glob(os.path.join(rundir, pattern)))
        if paths:
            return paths
    return []


def phase_rows(rundir):
    """Table-3 rows: per-tool phase means with their share of that tool's total."""
    rows = []

    malt = [malt_log_phases(p) for p in _logs(rundir, "logs/MALT/*.log")]
    malt = [m for m in malt if m["total"] > 0]
    if malt:
        total = stats.mean(m["total"] for m in malt)
        spec = [("load DB index from disk", "db_load", ""),
                ("alignment proper", "align", ""),
                ("write/copy SAM", "write_sam", ""),
                ("JVM start + taxonomy + LCA binning", "other",
                 "remainder; malt-run does not time these separately")]
        for label, key, note in spec:
            mean = stats.mean(m[key] for m in malt)
            rows.append({"tool": "malt", "phase": label, "mean_s": mean,
                         "share_pct": 100 * mean / total, "n": len(malt), "note": note})
        rows.append({"tool": "malt", "phase": "TOTAL", "mean_s": total,
                     "share_pct": 100.0, "n": len(malt),
                     "note": f"self-reported peak {stats.mean(m['peak_gb'] for m in malt):.1f} GB "
                             "(bounded by -Xmx, an upper bound not a requirement)"})

    strobe = [strobealign_log_phases(p) for p in _logs(
        rundir, "logs/ALIGNMENT/strobealign/*.log", "logs/ALIGNMENT/*.log")]
    strobe = [s for s in strobe if s["index"] or s["map"]]
    if strobe:
        wall = stats.mean(s["index"] + s["map"] for s in strobe)
        index_label = ("load prebuilt index" if all(s["prebuilt"] for s in strobe)
                       else "build index over reference FASTA")
        for label, key, note in [(index_label, "index", "wall"),
                                 ("mapping", "map", "wall, threaded")]:
            mean = stats.mean(s[key] for s in strobe)
            rows.append({"tool": "strobealign", "phase": label, "mean_s": mean,
                         "share_pct": 100 * mean / wall, "n": len(strobe), "note": note})
        for label, key in [("  of which finding hits", "hits_cpu"),
                           ("  of which chaining", "chain_cpu"),
                           ("  of which extend/pair seeds", "extend_cpu")]:
            rows.append({"tool": "strobealign", "phase": label,
                         "mean_s": stats.mean(s[key] for s in strobe), "share_pct": "",
                         "n": len(strobe), "note": "summed across threads, not wall"})
        rows.append({"tool": "strobealign", "phase": "TOTAL (index + mapping)",
                     "mean_s": wall, "share_pct": 100.0, "n": len(strobe),
                     "note": "the rule's wall time is larger: it also runs pigz + samtools sort"})

    # bowtie2 logs no timings, so it only appears through its ngsLCA step
    for tool, patterns in [("ngsLCA", ["logs/NGSLCA/*.log"]),
                           ("ngsLCA (strobealign)", ["logs/NGSLCA/strobealign/*.log"]),
                           ("ngsLCA (bowtie2)", ["logs/NGSLCA/bowtie2/*.log"]),
                           ("ngsLCA (bwa_aln)", ["logs/NGSLCA/bwa_aln/*.log"]),
                           ("ngsLCA (bwa_mem)", ["logs/NGSLCA/bwa_mem/*.log"])]:
        lca = [t for t in map(ngslca_log_walltime, _logs(rundir, *patterns))
               if t is not None]
        if lca:
            rows.append({"tool": tool, "phase": "TOTAL", "mean_s": stats.mean(lca),
                         "share_pct": 100.0, "n": len(lca), "note": "self-reported walltime"})
    return rows


# --------------------------------------------------------------------------- #
# Rendering
# --------------------------------------------------------------------------- #
def write_tsv(path, cols, rows, places=2):
    """Write dict rows as TSV, rounding floats -- these are timings measured once,
    so trailing precision is noise rather than information."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for row in rows:
            w.writerow([f"{v:.{places}f}" if isinstance(v, float) else v
                        for v in (row.get(c, "") for c in cols)])


def fmt_secs(seconds):
    return f"{seconds:.1f}s" + (f" ({seconds/60:.1f}m)" if seconds >= 120 else "")


def print_rule_table(title, summaries):
    if not summaries:
        return
    print(f"\n{title}")
    print(f"  {'rule':<28}{'n':>3}{'wall mean':>12}{'wall min':>10}{'wall max':>10}"
          f"{'peak RSS':>11}{'cpu time':>11}{'eff thr':>9}")
    print("  " + "-" * 92)
    for s in summaries:
        print(f"  {s['rule']:<28}{s['n']:>3}{s['wall_mean_s']:>11.1f}s"
              f"{s['wall_min_s']:>9.1f}s{s['wall_max_s']:>9.1f}s"
              f"{s['max_rss_max_gb']:>9.1f}GB{s['cpu_time_mean_s']:>10.1f}s"
              f"{s['eff_threads']:>9.1f}")


def print_profiler_table(by_profiler):
    if not by_profiler:
        return
    print("\nPER-SAMPLE COST BY PROFILER")
    print("  " + "-" * 92)
    totals = {}
    for profiler, entries in by_profiler.items():
        wall = sum(s["wall_mean_s"] for _stage, s in entries)
        rss = max(s["max_rss_max_gb"] for _stage, s in entries)
        totals[profiler] = (wall, rss)
        print(f"  {profiler}")
        for stage, s in entries:
            print(f"    {stage:<34}{s['wall_mean_s']:>9.1f}s{s['max_rss_max_gb']:>9.1f} GB"
                  f"   [{s['rule']}]")
        if len(entries) > 1:
            print(f"    {'TOTAL':<34}{wall:>9.1f}s{rss:>9.1f} GB")
        n = max(s["n"] for _stage, s in entries)
        print(f"    {f'x {n} samples':<34}{fmt_secs(wall * n):>16}")

    if len(totals) > 1:
        base, (wb, rb) = min(totals.items(), key=lambda kv: kv[1][0])
        for profiler, (wall, rss) in sorted(totals.items(), key=lambda kv: -kv[1][0]):
            if profiler != base and wb and rb:
                print(f"\n  {profiler} / {base}: {wall/wb:.2f}x wall, {rss/rb:.2f}x peak RSS",
                      end="")
        print()


def print_phase_table(rows):
    if not rows:
        return
    print("\nPHASE BREAKDOWN (from the tools' own logs, not from Snakemake)")
    print("  " + "-" * 92)
    tool = None
    for row in rows:
        if row["tool"] != tool:
            tool = row["tool"]
            print(f"  {tool}  (n={row['n']})")
        share = f"{row['share_pct']:>6.1f}%" if row["share_pct"] != "" else " " * 7
        print(f"    {row['phase']:<38}{row['mean_s']:>9.1f}s{share}"
              + (f"   {row['note']}" if row["note"] else ""))


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
MAP_COLOR = "#2a78d6"
INDEX_COLOR = "#eb6834"
LCA_COLOR = "#1baf7a"


def _bar_axes(plt, n_bars, ylabel, title):
    fig, ax = plt.subplots(figsize=(1.2 + 1.1 * n_bars, 3.6))
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", linewidth=0.4, color="0.88", zorder=0)
    return fig, ax


def _label_bars(ax, bars, values, fmt):
    for bar, v in zip(bars, values):
        ax.annotate(fmt(v), (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom", fontsize=8, color="0.25")


def _save(fig, outdir, name):
    fig.tight_layout()
    paths = [os.path.join(outdir, f"{name}.{ext}") for ext in ("pdf", "png")]
    for path in paths:
        fig.savefig(path, dpi=200)
    return paths


def plot_resources(summaries, outdir):
    """Three bar charts, one bar (group) per profiler: mean per-sample profiling
    wall time with the min-max over samples as a whisker, index build wall time,
    and peak RSS of profiling next to that of index construction.

    'Profiling' is all of a profiler's PROFILER_STAGES: alignment stacked under
    ngsLCA's taxonomic assignment, or for malt the one rule that does both. The
    whisker spans each sample's stages summed. A '*' on a profiler marks a rule
    read from --extra-benchmarks.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stages = {}
    for profiler, _stage, rule in PROFILER_STAGES:
        if rule in summaries:
            stages.setdefault(profiler, []).append(summaries[rule])
    index = {p: summaries[r] for p, r in INDEX_RULES.items() if r in summaries}
    written = []

    def tick(profiler, *recs):
        return profiler + ("*" if any(r["rule"].endswith("*") for r in recs) else "")

    if stages:
        names = list(stages)
        xs = range(len(names))
        n = max(s["n"] for p in names for s in stages[p])
        fig, ax = _bar_axes(plt, len(names), "wall time per sample (min)",
                            f"Profiling runtime, mean over {n} samples")
        bottom = [0.0] * len(names)
        for i, (color, label) in enumerate(((MAP_COLOR, "alignment"),
                                            (LCA_COLOR, "taxonomic assignment"))):
            height = [stages[p][i]["wall_mean_s"] / 60 if i < len(stages[p]) else 0
                      for p in names]
            if any(height):
                # white edge = 2px surface gap between stacked segments
                ax.bar(xs, height, bottom=bottom, width=0.6, color=color, label=label,
                       edgecolor="white", linewidth=1, zorder=2)
                bottom = [b + h for b, h in zip(bottom, height)]
        lo, hi = [], []
        for p, top in zip(names, bottom):
            jobs = set.intersection(*(set(s["wall_by_job"]) for s in stages[p]))
            totals = [sum(s["wall_by_job"][j] for s in stages[p]) / 60 for j in jobs]
            lo.append(max(0.0, top - min(totals)) if totals else 0.0)
            hi.append(max(0.0, max(totals) - top) if totals else 0.0)
        ax.errorbar(xs, bottom, yerr=[lo, hi], fmt="none", elinewidth=1, capsize=3,
                    ecolor="0.35", zorder=3)
        for x, top, up in zip(xs, bottom, hi):    # labels clear the whisker
            ax.annotate(f"{top:.1f}", (x, top + up), xytext=(0, 3),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=8, color="0.25")
        ax.set_xticks(xs, [tick(p, *stages[p]) for p in names])
        ax.margins(y=0.12)
        ax.legend(frameon=False, fontsize=8, loc="upper right")
        written += _save(fig, outdir, "bench_profiling_runtime")
        plt.close(fig)

    if index:
        names = list(index)
        hours = [index[p]["wall_mean_s"] / 3600 for p in names]
        fig, ax = _bar_axes(plt, len(names), "wall time (h)", "Index construction runtime")
        bars = ax.bar(range(len(names)), hours, width=0.6, color=INDEX_COLOR, zorder=2)
        _label_bars(ax, bars, hours,
                    lambda v: f"{v:.2f}" if v < 1 else f"{v:.1f}")
        ax.set_xticks(range(len(names)), [tick(p, index[p]) for p in names])
        ax.margins(y=0.12)
        written += _save(fig, outdir, "bench_index_runtime")
        plt.close(fig)

    profiling = {p: max(ss, key=lambda s: s["max_rss_max_gb"]) for p, ss in stages.items()}
    names = [p for p in INDEX_RULES if p in profiling or p in index]
    if names:
        width = 0.38
        fig, ax = _bar_axes(plt, len(names), "peak RSS (GB)", "Peak memory")
        for offset, recs, color, label in (
                (-width / 2, profiling, MAP_COLOR, "profiling"),
                (width / 2, index, INDEX_COLOR, "index construction")):
            xs = [i + offset for i, p in enumerate(names) if p in recs]
            gb = [recs[p]["max_rss_max_gb"] for p in names if p in recs]
            # 2px surface gap between the paired bars
            bars = ax.bar(xs, gb, width=width, color=color, label=label,
                          edgecolor="white", linewidth=1, zorder=2)
            _label_bars(ax, bars, gb, lambda v: f"{v:.1f}")
        ax.set_xticks(range(len(names)),
                      [tick(p, *(r[p] for r in (profiling, index) if p in r)) for p in names])
        ax.margins(y=0.12)
        ax.legend(frameon=False, fontsize=8, loc="upper right")
        written += _save(fig, outdir, "bench_peak_memory")
        plt.close(fig)
    return written


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rundir", default=".",
                    help="aMeta run directory holding benchmarks/ and logs/ (default: .)")
    ap.add_argument("--outdir", default="",
                    help="also write bench_rules.tsv, bench_phases.tsv and the runtime/memory "
                         "plots here")
    ap.add_argument("--rules", default="",
                    help="comma-separated extra rule names to include in the "
                         "per-rule table (default: the profiler, setup and "
                         "shared rules this script knows about)")
    ap.add_argument("--extra-benchmarks", default="",
                    help="benchmarks dir to read rules from when the run has none "
                         "for them, e.g. KrakenUniq reused from an earlier run")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    bench_dir = os.path.join(args.rundir, "benchmarks")
    if not os.path.isdir(bench_dir):
        sys.exit(f"No benchmarks/ under {args.rundir} -- pass --rundir <aMeta run dir>")

    profiler_rules = [rule for _p, _s, rule in PROFILER_STAGES]
    extra = [r for r in args.rules.split(",") if r]
    wanted = profiler_rules + list(SETUP_RULES) + list(SHARED_RULES) + extra

    bench_dirs = [bench_dir] + ([args.extra_benchmarks] if args.extra_benchmarks else [])
    summaries, missing = {}, []
    for rule in wanted:
        s = summarise_rule(bench_dirs, rule)
        if s:
            summaries[rule] = s
        else:
            missing.append(rule)

    print(f"Benchmarks under {os.path.abspath(bench_dir)}")
    if args.extra_benchmarks:
        print(f"  * from {os.path.abspath(args.extra_benchmarks)}")
    if missing:
        print(f"  (no benchmark files for: {', '.join(missing)})")

    print_rule_table("PER SAMPLE", [summaries[r] for r in profiler_rules if r in summaries])
    print_rule_table("ONE-OFF SETUP", [summaries[r] for r in SETUP_RULES if r in summaries])
    print_rule_table("SHARED UPSTREAM (paid by every profiler)",
                     [summaries[r] for r in SHARED_RULES if r in summaries])
    if extra:
        print_rule_table("ALSO REQUESTED", [summaries[r] for r in extra if r in summaries])

    by_profiler = {}
    for profiler, stage, rule in PROFILER_STAGES:
        if rule in summaries:
            by_profiler.setdefault(profiler, []).append((stage, summaries[rule]))
    print_profiler_table(by_profiler)

    phases = phase_rows(args.rundir)
    print_phase_table(phases)

    print("\nRead these with the caveats in the module docstring: wall times only "
          "compare\nwhen the jobs ran serially, and a JVM tool's peak RSS is bounded "
          "by its -Xmx heap.")

    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)
        rules_tsv = os.path.join(args.outdir, "bench_rules.tsv")
        phases_tsv = os.path.join(args.outdir, "bench_phases.tsv")
        write_tsv(rules_tsv, RULE_COLS, [summaries[r] for r in wanted if r in summaries])
        write_tsv(phases_tsv, PHASE_COLS, phases)
        print(f"\n{rules_tsv}\n{phases_tsv}")
        for path in plot_resources(summaries, args.outdir):
            print(path)


if __name__ == "__main__":
    main()
