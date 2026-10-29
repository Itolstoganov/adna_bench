#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared" / "scripts"))
import tool_style  # noqa: E402

METRICS = [
    ("% Mapped endogenous reads", "mapped_endogenous"),
    ("% Accuracy (endogenous only)", "accuracy"),
    ("% Mapped deaminated endogenous reads", "mapped_endogenous_deam"),
    ("% Accuracy (deaminated endogenous only)", "accuracy_deam"),
    ("% Mapped bacterial reads", "mapped_bacterial"),
    ("% Mapped contaminated reads", "mapped_contaminated"),
]

# Optional secondary-alignment accuracy. Detected by column presence in the table.
SECONDARY_METRICS = [
    ("% Accuracy with secondary (endogenous only)", "accuracy_secondary"),
    ("% Accuracy with secondary (deaminated endogenous only)", "accuracy_secondary_deam"),
]


def with_secondary_if_present(base_metrics, df):
    """Append the secondary metrics whose columns exist in `df` (enabled runs only)."""
    return base_metrics + [m for m in SECONDARY_METRICS if m[0] in df.columns]

# Resource metrics are k-metrics only (line vs k); column, file slug, plot title.
RESOURCE_METRICS = [
    ("Mapping Time (s)", "k_runtime", "Mapping Time"),
    ("Peak Memory (GB)", "k_memory", "Peak Memory"),
]

COMBINED_K_METRICS = [
    ("% Accuracy (endogenous only)", "accuracy"),
    ("% Accuracy (deaminated endogenous only)", "accuracy_deam"),
    ("% Mapped endogenous reads", "mapped_endogenous"),
    ("% Mapped deaminated endogenous reads", "mapped_endogenous_deam"),
]


# --- filters -----------------------------------------------------------------
# Filters are callables df -> df. They are column-tolerant so the same filter
# can be applied to the accuracy table and the resource table (which lacks a
# "Score Threshold" column).

def score_eq(threshold):
    """Default filter: keep rows whose Score Threshold equals `threshold`."""
    def f(df):
        if "Score Threshold" not in df.columns:
            return df
        return df[df["Score Threshold"] == threshold]
    return f


def k_eq(k):
    """Keep rows for a single k value, plus k-less tools (k is NaN)."""
    def f(df):
        return df[(df["k"] == k) | (df["k"].isna())]
    return f


# --- helpers -----------------------------------------------------------------

# Strobealign version names in config order (via --strobealign); they pin each
# version to the same color and marker in every plot, whichever subset of tools
# a plot shows. Aligner styles are fixed in shared/scripts/tool_style.py.
STROBEALIGN_VERSIONS = []


def _palette(names):
    """Tool -> color, shared with the metagenomic benchmark."""
    return {n: c for n, (c, _m) in tool_style.styles(names, STROBEALIGN_VERSIONS).items()}


def _markers(names):
    """Tool -> marker, shared with the metagenomic benchmark."""
    return {n: m for n, (_c, m) in tool_style.styles(names, STROBEALIGN_VERSIONS).items()}


def read_table(path):
    """An accuracy or resource table with aligner names in display form.

    The tables keep the raw program names ("bwaaln", "bwa mem (k=17)") so they
    do not depend on presentation; renaming happens here, before any plotting.
    """
    df = pd.read_csv(path, sep="\t")
    for col in ("Tool", "Base Name"):
        if col in df.columns:
            df[col] = df[col].map(lambda n: n if pd.isna(n) else tool_style.display(n))
    return df


def _split_baseline(df, x, baseline):
    """Split rows that have a value on the x axis from k-less baseline rows.

    Returns (varying_df, baseline_df). When `baseline` is False or x has no NaN
    values, baseline_df is empty and varying_df contains everything.
    """
    if baseline and df[x].isna().any():
        varying = df[df[x].notna()].copy()
        base = df[df[x].isna()].copy()
    else:
        varying = df.copy()
        base = df.iloc[0:0].copy()
    if not varying.empty and x == "k":
        varying["k"] = varying["k"].astype(int)
    return varying, base


# --- unified plotting core ---------------------------------------------------

def plot_xy(df, x, y, *, output_path=None, pdf=None, filter=None,
            hue: str | None = "Tool", baseline_hue="Tool",
            baseline=True, row=None, col=None, title=None):
    """Plot line `y` against column `x`, grouped by `hue`, after `filter`.

    - `filter`: callable df -> df (default: identity).
    - `baseline`: draw rows whose `x` is NaN (k-less tools on a k axis) as
      horizontal dashed reference lines (single-axes) or dashed broadcast lines
      (faceted), colored/labeled by `baseline_hue`.
    - `row`/`col`: facet columns; when set, render with sns.relplot and write to
      `pdf` (a PdfPages) if given, else to `output_path`.
    - exactly one of `output_path` / `pdf` is used as the output target.
    """
    if filter is not None:
        df = filter(df)
    if df.empty:
        return

    if row is not None or col is not None:
        _plot_faceted(df, x, y, output_path=output_path, pdf=pdf, hue=hue,
                      baseline=baseline, row=row, col=col, title=title)
        return

    varying, base = _split_baseline(df, x, baseline)
    if varying.empty and base.empty:
        return

    names = list(varying[hue].unique()) if hue else []
    names += list(base[baseline_hue].unique())
    palette = _palette(names) if names else None
    markers = _markers(names) if names else None

    plt.figure(figsize=(10, 6))
    if not varying.empty:
        if hue:
            sns.lineplot(data=varying, x=x, y=y, hue=hue, palette=palette,
                         style=hue, markers=markers, dashes=False)
        else:
            sns.lineplot(data=varying, x=x, y=y, marker="o")
        plt.xticks(sorted(varying[x].unique()))

    for _, r in base.iterrows():
        name = r[baseline_hue]
        plt.axhline(y=r[y], linestyle="--", color=palette[name], label=name)

    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(title if title is not None else y)
    if hue or not base.empty:
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def _plot_faceted(df, x, y, *, output_path, pdf, hue, baseline, row, col, title):
    """Faceted (row x col) line plot; k-less tools broadcast as dashed lines."""
    is_k_aware = df[x].notna()
    df = df.copy()
    df["k_aware"] = is_k_aware

    k_aware_df = df[df["k_aware"]].copy()
    if k_aware_df.empty:
        return
    k_aware_df[x] = k_aware_df[x].astype(int)
    x_values = sorted(k_aware_df[x].unique())

    base_df = df[~df["k_aware"]].copy() if baseline else df.iloc[0:0].copy()
    if not base_df.empty and x_values:
        broadcast = base_df.loc[base_df.index.repeat(len(x_values))].copy()
        broadcast[x] = x_values * len(base_df)
        plot_df = pd.concat([k_aware_df, broadcast], ignore_index=True)
    else:
        plot_df = k_aware_df

    display_order = []
    for name in list(k_aware_df[hue].unique()) + list(base_df[hue].unique()):
        if name not in display_order:
            display_order.append(name)
    palette = _palette(display_order)
    markers = _markers(display_order)
    base_displays = set(base_df[hue].unique())
    dashes = {d: (4, 2) if d in base_displays else "" for d in display_order}

    facet_kwargs = {}
    if row is not None:
        facet_kwargs["row"] = row
        facet_kwargs["row_order"] = sorted(plot_df[row].unique())
    if col is not None:
        facet_kwargs["col"] = col
        facet_kwargs["col_order"] = sorted(plot_df[col].unique())

    g = sns.relplot(
        data=plot_df,
        x=x, y=y,
        hue=hue, hue_order=display_order, palette=palette,
        style=hue, style_order=display_order, markers=markers, dashes=dashes,
        kind="line",
        facet_kws={"sharey": False, "legend_out": True},
        legend=False,
        **facet_kwargs,
    )
    g.set_axis_labels(x, y)
    g.set_titles(row_template="{row_name}", col_template="{col_name}")
    for ax in g.axes.flat:
        ax.set_xticks(x_values)

    handles = []
    for display in display_order:
        k_aware = display not in base_displays
        handles.append(Line2D(
            [], [],
            color=palette[display],
            linestyle="-" if k_aware else "--",
            marker=markers[display],
        ))
    g.figure.legend(
        handles, display_order, title="Tool",
        loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False,
    )

    g.figure.suptitle(title if title is not None else y)
    g.tight_layout()
    if pdf is not None:
        pdf.savefig(g.figure, bbox_inches="tight")
    else:
        g.figure.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(g.figure)


def _build_combined_frame(triples, score_threshold=None):
    """Concatenate per-dataset tables tagged with facet row/col labels.

    `triples` is a list of (path, row_label, col_label). When `score_threshold`
    is given and the table has a "Score Threshold" column (accuracy tables), rows
    are filtered to that threshold; resource tables (no such column) pass through.
    """
    frames = []
    for path, row_label, col_label in triples:
        df = read_table(path)
        if score_threshold is not None and "Score Threshold" in df.columns:
            df = df[df["Score Threshold"] == score_threshold]
        df = df.copy()
        df["dataset_row"] = row_label
        df["dataset_col"] = col_label
        frames.append(df)
    if not frames:
        return None
    df = pd.concat(frames, ignore_index=True)
    is_k_aware = df["k"].astype(str).str.strip() != "nan"
    df["Display"] = df["Base Name"].where(is_k_aware, df["Tool"])
    return df


def plot_combined_k_pdf(dataset_inputs, output_path, score_threshold):
    # dataset_inputs: list of (NAME, ACC_PATH, RES_PATH, ROW, COL)
    acc_df = _build_combined_frame(
        [(acc, row, col) for _, acc, _res, row, col in dataset_inputs],
        score_threshold,
    )
    res_df = _build_combined_frame(
        [(res, row, col) for _, _acc, res, row, col in dataset_inputs],
    )
    combined_metrics = COMBINED_K_METRICS
    if acc_df is not None:
        combined_metrics = with_secondary_if_present(combined_metrics, acc_df)

    with PdfPages(output_path) as pdf:
        if acc_df is not None and not acc_df.empty:
            for metric_col, _ in combined_metrics:
                plot_xy(acc_df, x="k", y=metric_col, pdf=pdf,
                        hue="Display", baseline=True,
                        row="dataset_row", col="dataset_col", title=metric_col)
        if res_df is not None and not res_df.empty:
            for y_col, _file_name, title in RESOURCE_METRICS:
                plot_xy(res_df, x="k", y=y_col, pdf=pdf,
                        hue="Display", baseline=True,
                        row="dataset_row", col="dataset_col", title=title)


# --- summary plots (datasets on x axis) --------------------------------------

def _build_short_names(triples):
    """Build shortened dataset labels `genome|read_len|variation`.

    `triples` is a list of (genome, read_length, variation). The read length and
    variation components are omitted when every dataset shares the same value, so
    e.g. with a single read length the label becomes `CHM13|var3`.
    """
    read_lens = {rl for _g, rl, _v in triples}
    variations = {v for _g, _rl, v in triples}
    include_rl = len(read_lens) > 1
    include_var = len(variations) > 1
    names = []
    for genome, read_len, variation in triples:
        parts = [str(genome)]
        if include_rl:
            parts.append(str(read_len))
        if include_var:
            parts.append(str(variation))
        names.append("|".join(parts))
    return names


def _build_summary_frame(pairs, *, k=None, score_threshold=None):
    """Concatenate per-dataset tables tagged with a short dataset label.

    `pairs` is a list of (path, short_name). Rows are filtered to a single k
    (keeping k-less tools, whose k is NaN) and, for accuracy tables, to a single
    score threshold.
    """
    frames = []
    for path, short_name in pairs:
        df = read_table(path)
        if score_threshold is not None and "Score Threshold" in df.columns:
            df = df[df["Score Threshold"] == score_threshold]
        if k is not None:
            df = df[(df["k"] == k) | (df["k"].isna())]
        df = df.copy()
        df["dataset"] = short_name
        frames.append(df)
    if not frames:
        return None
    df = pd.concat(frames, ignore_index=True)
    is_k_aware = df["k"].notna()
    df["Display"] = df["Base Name"].where(is_k_aware, df["Tool"])
    return df


# Tool treated as the reference for delta (difference) plots.
SUMMARY_DELTA_BASELINE = tool_style.display("bwaaln")

# Friendlier metric labels for summary plot titles / axes. Keys are the table
# column names; everything not listed is shown verbatim.
SUMMARY_METRIC_LABELS = {
    "% Accuracy (endogenous only)": "% Accuracy",
    "% Accuracy (deaminated endogenous only)": "% Accuracy (deaminated reads only)",
}


def _summary_label(col):
    return SUMMARY_METRIC_LABELS.get(col, col)


def _apply_baseline(df, y, baseline, op):
    """Rebase `y` against `baseline`'s per-dataset value using `op`.

    `op` is "sub" (difference) or "div" (ratio). Datasets where the baseline is
    absent get NaN and are dropped downstream.
    """
    base = (df[df["Display"] == baseline][["dataset", y]]
            .rename(columns={y: "__base"}))
    merged = df.merge(base, on="dataset", how="left").copy()
    if op == "div":
        merged[y] = merged[y] / merged["__base"]
    else:
        merged[y] = merged[y] - merged["__base"]
    return merged.drop(columns="__base")


def plot_summary_metric(df, y, dataset_order, *, pdf=None, png_path=None,
                        title=None, delta_baseline=None, ratio_baseline=None,
                        log_y=False, x="dataset", xlabel="Dataset"):
    """Point/line plot of metric `y` with datasets (or column `x`) on the x axis,
    one line per tool.

    When `delta_baseline` is set, values are plotted as the difference to that
    tool (a dashed line marks zero); when `ratio_baseline` is set, as the ratio to
    that tool (a dashed line marks one). `log_y` puts the y axis on a log scale.
    The figure is written to `pdf` and/or `png_path` (both may be given).
    """
    if y not in df.columns:
        return
    if delta_baseline is not None:
        df = _apply_baseline(df, y, delta_baseline, "sub")
    elif ratio_baseline is not None:
        df = _apply_baseline(df, y, ratio_baseline, "div")
    sub = df[df[y].notna()]
    if sub.empty:
        return
    names = list(dict.fromkeys(sub["Display"]))
    palette = _palette(names)
    markers = _markers(names)
    plt.figure(figsize=(max(10, len(dataset_order) * 1.2), 6))
    if delta_baseline is not None:
        plt.axhline(0, color="gray", linewidth=0.8, linestyle="--")
    elif ratio_baseline is not None:
        plt.axhline(1, color="gray", linewidth=0.8, linestyle="--")
    sns.pointplot(data=sub, x=x, y=y, hue="Display",
                  order=dataset_order, hue_order=names, palette=palette,
                  markers=[markers[n] for n in names], errorbar=None)
    label = _summary_label(y)
    if delta_baseline is not None:
        ylabel = f"Δ {label} (vs {delta_baseline})"
    elif ratio_baseline is not None:
        ylabel = f"{label} (ratio to {ratio_baseline})"
    else:
        ylabel = label
    if log_y:
        plt.yscale("log")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title if title is not None else y)
    plt.xticks(rotation=45, ha="right")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    if pdf is not None:
        pdf.savefig(bbox_inches="tight")
    if png_path is not None:
        plt.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_summary_pdf(summary_inputs, output_path, k, score_threshold,
                     png_dir=None, num_reads=None):
    # summary_inputs: list of (NAME, ACC_PATH, RES_PATH, GENOME, READ_LEN, VARIATION, GT_PATH)
    short_names = _build_short_names(
        [(genome, read_len, var) for _n, _a, _r, genome, read_len, var, _gt in summary_inputs]
    )
    dataset_order = list(dict.fromkeys(short_names))
    genome_by_dataset = {
        short: genome
        for (_n, _a, _r, genome, _rl, _v, _gt), short in zip(summary_inputs, short_names)
    }
    read_len_by_dataset = {
        short: float(read_len)
        for (_n, _a, _r, _g, read_len, _v, _gt), short in zip(summary_inputs, short_names)
    }
    variation_by_dataset = {
        short: _variation_label(var)
        for (_n, _a, _r, _g, _rl, var, _gt), short in zip(summary_inputs, short_names)
    }
    deam_fraction_by_dataset = {
        short: count_deaminated_fraction(gt)
        for (_n, _a, _r, _g, _rl, _v, gt), short in zip(summary_inputs, short_names)
    }
    acc_df = _build_summary_frame(
        [(acc, short) for (_n, acc, _r, *_), short in zip(summary_inputs, short_names)],
        k=k, score_threshold=score_threshold,
    )
    res_df = _build_summary_frame(
        [(res, short) for (_n, _a, res, *_), short in zip(summary_inputs, short_names)],
        k=k,
    )
    metrics = METRICS
    if acc_df is not None:
        metrics = with_secondary_if_present(metrics, acc_df)

    if png_dir is not None:
        png_dir.mkdir(parents=True, exist_ok=True)

    def png_path(slug):
        return png_dir / f"{slug}.png" if png_dir is not None else None

    variation_order = sorted(set(variation_by_dataset.values()), key=_variation_sort_key)

    def with_reference_columns(df):
        df = df.copy()
        df["reference"] = df["dataset"].map(genome_by_dataset)
        df["variation"] = df["dataset"].map(variation_by_dataset)
        return df

    with PdfPages(output_path) as pdf:
        if acc_df is not None and not acc_df.empty:
            # absolute values
            for metric_col, slug in metrics:
                plot_summary_metric(acc_df, metric_col, dataset_order, pdf=pdf,
                                    png_path=png_path(slug),
                                    title=_summary_label(metric_col))
            # delta vs the default tool
            for metric_col, slug in metrics:
                plot_summary_metric(
                    acc_df, metric_col, dataset_order, pdf=pdf,
                    png_path=png_path(f"{slug}_delta"),
                    title=f"{_summary_label(metric_col)} (Δ vs {SUMMARY_DELTA_BASELINE})",
                    delta_baseline=SUMMARY_DELTA_BASELINE)
        if res_df is not None and not res_df.empty:
            # runtime normalised per processed megabase: total reads (config) x
            # read length (per dataset). Falls back to total time if num_reads is unset.
            runtime_col = "Mapping Time (s)"
            runtime_title = "Mapping Time"
            if num_reads:
                res_df = res_df.copy()
                read_len = res_df["dataset"].map(read_len_by_dataset)
                res_df["Runtime per Mbp (s)"] = (
                    res_df["Mapping Time (s)"] / (num_reads * read_len / 1e6)
                )
                runtime_col = "Runtime per Mbp (s)"
                runtime_title = "Runtime per Mbp"
            plot_reference_lines(with_reference_columns(res_df), runtime_col,
                                 variation_order, pdf=pdf,
                                 png_path=png_path("k_runtime"), title=runtime_title,
                                 sharey=False)
            # absolute peak memory
            plot_summary_metric(res_df, "Peak Memory (GB)", dataset_order, pdf=pdf,
                                png_path=png_path("k_memory"), title="Peak Memory")
            # peak memory relative to the default aligner (log scale)
            plot_summary_metric(
                res_df, "Peak Memory (GB)", dataset_order, pdf=pdf,
                png_path=png_path("k_memory_delta"),
                title=f"Peak Memory (relative to {SUMMARY_DELTA_BASELINE})",
                ratio_baseline=SUMMARY_DELTA_BASELINE, log_y=True)
            plot_index_metrics(res_df, genome_by_dataset, pdf=pdf, png_path=png_path)
        if acc_df is not None and not acc_df.empty:
            acc_df = with_reference_columns(acc_df)
            plot_outcome_bars(acc_df, "% Accuracy (endogenous only)",
                              "% Mapped endogenous reads", variation_order,
                              pdf=pdf, png_path=png_path("outcome"),
                              title="Mapping accuracy")
            plot_outcome_bars(acc_df, "% Accuracy (deaminated endogenous only)",
                              "% Mapped deaminated endogenous reads", variation_order,
                              pdf=pdf, png_path=png_path("outcome_deam"),
                              title="Mapping accuracy (deaminated reads)")
            acc_df["deam_fraction"] = acc_df["dataset"].map(deam_fraction_by_dataset)
            plot_deam_accuracy_diff(acc_df, variation_order, pdf=pdf,
                                    png_path=png_path("accuracy_deam_diff"))


def _variation_label(variation):
    """`var0p1` -> `0.1%`; unrecognised values are returned unchanged."""
    try:
        return f"{float(str(variation).removeprefix('var').replace('p', '.')):g}%"
    except ValueError:
        return str(variation)


def _variation_sort_key(label):
    try:
        return (0, float(label.removesuffix("%")), label)
    except ValueError:
        return (1, 0.0, label)


def count_deaminated_fraction(gt_path):
    """Fraction of endogenous reads in the ground truth that carry deamination."""
    endo = deam = 0
    with open(gt_path) as fh:
        for line in fh:
            name = line.split("\t", 4)[3]
            origin = name.split(":", 4)[4].split("_DEAM:", 1)[0].lstrip("0123456789")
            if origin.startswith("e"):
                endo += 1
                deam += "_DEAM:" in name
    return deam / endo if endo else float("nan")


def _reference_axes(references, *, sharey, height=5.5):
    fig, axes = plt.subplots(1, len(references), figsize=(6 * len(references), height),
                             sharey=sharey, squeeze=False)
    return fig, axes[0]


def _save_summary_figure(fig, pdf, png_path):
    fig.tight_layout()
    if pdf is not None:
        pdf.savefig(fig, bbox_inches="tight")
    if png_path is not None:
        fig.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_outcome_bars(df, correct_col, mapped_col, variation_order, *, pdf=None,
                      png_path=None, title=None):
    """Stacked correct / incorrect / unmapped bars in a single panel: variation
    levels grouped by reference on the x axis, one bar per tool within each level."""
    if correct_col not in df.columns or mapped_col not in df.columns:
        return
    df = df[df[correct_col].notna() & df[mapped_col].notna()]
    if df.empty:
        return
    references = list(dict.fromkeys(df["reference"]))
    tools = list(dict.fromkeys(df["Display"]))
    palette = _palette(tools)
    bar_width = 0.8 / len(tools)

    groups = []
    for reference in references:
        present = set(df.loc[df["reference"] == reference, "variation"])
        groups.append((reference, [v for v in variation_order if v in present]))
    reference_gap = 0.6
    centers, labels, spans = [], [], []
    x = 0.0
    for reference, levels in groups:
        first = x
        for level in levels:
            centers.append(x)
            labels.append((reference, level))
            x += 1
        spans.append((reference, first, x - 1))
        x += reference_gap

    fig, ax = plt.subplots(figsize=(max(12, 0.42 * len(centers) * len(tools)), 4))
    for i, tool in enumerate(tools):
        rows = df[df["Display"] == tool].set_index(["reference", "variation"])
        xs, correct, incorrect, unmapped = [], [], [], []
        for center, key in zip(centers, labels):
            if key not in rows.index:
                continue
            r = rows.loc[key]
            xs.append(center - 0.4 + bar_width * (i + 0.5))
            correct.append(r[correct_col])
            incorrect.append(r[mapped_col] - r[correct_col])
            unmapped.append(100 - r[mapped_col])
        color = palette[tool]
        ax.bar(xs, correct, bar_width * 0.9, color=color, edgecolor="white", linewidth=0.5)
        ax.bar(xs, incorrect, bar_width * 0.9, bottom=correct, color=color, alpha=0.4,
               hatch="///", edgecolor="white", linewidth=0.5)
        ax.bar(xs, unmapped, bar_width * 0.9,
               bottom=[c + w for c, w in zip(correct, incorrect)],
               color="#e5e7eb", edgecolor="white", linewidth=0.5)
        for bx, c, w in zip(xs, correct, incorrect):
            ax.text(bx, 101, f"{c:.1f}", ha="center", va="bottom", fontsize=6.5)
            ax.text(bx, 107, f"{c + w:.1f}", ha="center", va="bottom", fontsize=6.5,
                    color="#6b7280")
    ax.set_xticks(centers)
    ax.set_xticklabels([level for _ref, level in labels])
    ax.set_xlim(centers[0] - 0.5, centers[-1] + 0.5)
    for (_ref, _first, last), (_next, next_first, _l) in zip(spans, spans[1:]):
        ax.axvline((last + next_first) / 2, color="#9ca3af", linewidth=0.8)
    for reference, first, last in spans:
        ax.text((first + last) / 2, -0.12, reference, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=11, fontweight="bold")
    ax.set_ylim(0, 115)
    ax.set_yticks(range(0, 101, 20))
    ax.set_ylabel("Reads (%)")
    ax.spines[["top", "right"]].set_visible(False)

    tool_handles = [Patch(facecolor=palette[t]) for t in tools]
    outcome_handles = [
        Patch(facecolor="#6b7280"),
        Patch(facecolor="#6b7280", alpha=0.4, hatch="///", edgecolor="white"),
        Patch(facecolor="#e5e7eb"),
    ]
    tool_legend = fig.legend(tool_handles, tools, title="Tool", loc="center left",
                             bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.add_artist(tool_legend)
    fig.legend(outcome_handles, ["Correct", "Incorrect", "Unmapped"],
               title="Mapping outcome", loc="upper center", bbox_to_anchor=(0.5, -0.02),
               ncol=3, frameon=False)
    fig.suptitle(title if title is not None else correct_col)
    _save_summary_figure(fig, pdf, png_path)


def plot_deam_accuracy_diff(df, variation_order, *, pdf=None, png_path=None):
    """Accuracy on deaminated minus accuracy on non-deaminated endogenous reads,
    one panel per reference with variation levels on the x axis.

    Non-deaminated accuracy is recovered from the overall and deaminated
    accuracies using the fraction `f` of deaminated endogenous reads:
    acc_nondeam = (acc_all - f * acc_deam) / (1 - f).
    """
    all_col = "% Accuracy (endogenous only)"
    deam_col = "% Accuracy (deaminated endogenous only)"
    if all_col not in df.columns or deam_col not in df.columns:
        return
    df = df.copy()
    f = df["deam_fraction"]
    nondeam = (df[all_col] - f * df[deam_col]) / (1 - f)
    y = "Δ % Accuracy (deaminated − non-deaminated)"
    df[y] = df[deam_col] - nondeam
    plot_reference_lines(df, y, variation_order, pdf=pdf, png_path=png_path,
                         title="Accuracy difference: deaminated vs non-deaminated reads",
                         zero_line=True)


def plot_reference_lines(df, y, variation_order, *, pdf=None, png_path=None,
                         title=None, zero_line=False, sharey=True):
    """Point/line plot of `y`, one panel per reference, variation levels on the
    x axis and one line per tool."""
    if y not in df.columns:
        return
    df = df[df[y].notna()]
    if df.empty:
        return
    references = list(dict.fromkeys(df["reference"]))
    tools = list(dict.fromkeys(df["Display"]))
    palette = _palette(tools)
    markers = _markers(tools)

    fig, axes = _reference_axes(references, sharey=sharey)
    for ax, reference in zip(axes, references):
        ref_df = df[df["reference"] == reference]
        levels = [v for v in variation_order if v in set(ref_df["variation"])]
        if zero_line:
            ax.axhline(0, color="gray", linewidth=0.8, linestyle="--")
        sns.pointplot(data=ref_df, x="variation", y=y, hue="Display",
                      order=levels, hue_order=tools, palette=palette,
                      markers=[markers[t] for t in tools], errorbar=None, ax=ax,
                      legend=False)
        ax.set_title(reference)
        ax.set_xlabel("Variation")
        ax.set_ylabel(y)
        ax.yaxis.set_tick_params(labelleft=True)
    handles = [Line2D([], [], color=palette[t], marker=markers[t]) for t in tools]
    fig.legend(handles, tools, title="Tool", loc="center left",
               bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.suptitle(title if title is not None else y)
    _save_summary_figure(fig, pdf, png_path)


# Index construction metrics: column, file slug, plot title.
INDEX_METRICS = [
    ("Index Time (s)", "index_runtime", "Index Construction Time"),
    ("Index Peak Memory (GB)", "index_memory", "Index Construction Peak Memory"),
]


def plot_index_metrics(res_df, genome_by_dataset, *, pdf=None, png_path=None):
    """Index time/memory with references on the x axis.

    The index depends only on the reference, so datasets sharing a reference
    (differing in variation or read length) are collapsed to their median.
    """
    columns = [col for col, _slug, _title in INDEX_METRICS if col in res_df.columns]
    if not columns:
        return
    df = res_df.copy()
    df["reference"] = df["dataset"].map(genome_by_dataset)
    df = df.groupby(["reference", "Display"], as_index=False, sort=False)[columns].median()
    reference_order = list(dict.fromkeys(genome_by_dataset.values()))
    for col, slug, title in INDEX_METRICS:
        if col in columns:
            plot_summary_metric(df, col, reference_order, pdf=pdf,
                                png_path=png_path(slug), title=title,
                                x="reference", xlabel="Reference")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path)
    parser.add_argument("--resources", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--k-plot", action="store_true")
    parser.add_argument("--score-threshold", type=int, default=60)
    parser.add_argument("--k", type=int)
    parser.add_argument("--combined-pdf", type=Path,
                        help="Write a single combined PDF with k-plots for all "
                             "datasets and selected metrics, then exit.")
    parser.add_argument("--dataset", nargs=5, action="append", default=[],
                        metavar=("NAME", "ACC_PATH", "RES_PATH", "ROW", "COL"),
                        help="Dataset name, accuracy table path, resource table "
                             "path, row label, column label (repeatable, used "
                             "with --combined-pdf).")
    parser.add_argument("--summary-pdf", type=Path,
                        help="Write a single PDF with one summary plot per metric, "
                             "datasets on the x axis (at --main-k), then exit.")
    parser.add_argument("--summary-png-dir", type=Path,
                        help="Also write each summary plot as a PNG into this "
                             "directory (used with --summary-pdf).")
    parser.add_argument("--main-k", type=int,
                        help="k value selected for the summary plots.")
    parser.add_argument("--num-reads", type=int,
                        help="Total reads per dataset (config); enables runtime-"
                             "per-base in the summary plots.")
    parser.add_argument("--summary-dataset", nargs=7, action="append", default=[],
                        metavar=("NAME", "ACC_PATH", "RES_PATH", "GENOME",
                                 "READ_LEN", "VARIATION", "GT_PATH"),
                        help="Dataset name, accuracy table path, resource table "
                             "path, genome, read length, variation, ground truth "
                             "BED (repeatable, used with --summary-pdf).")
    parser.add_argument("--strobealign", action="append", default=[],
                        dest="strobealign_versions",
                        help="Strobealign version name, in config order "
                             "(repeatable); pins its color and marker across plots.")
    args = parser.parse_args()
    STROBEALIGN_VERSIONS.extend(args.strobealign_versions)

    if args.combined_pdf:
        if not args.dataset:
            parser.error("--combined-pdf requires at least one --dataset entry")
        args.combined_pdf.parent.mkdir(parents=True, exist_ok=True)
        plot_combined_k_pdf(args.dataset, args.combined_pdf, args.score_threshold)
        return

    if args.summary_pdf:
        if not args.summary_dataset:
            parser.error("--summary-pdf requires at least one --summary-dataset entry")
        args.summary_pdf.parent.mkdir(parents=True, exist_ok=True)
        plot_summary_pdf(args.summary_dataset, args.summary_pdf,
                         args.main_k, args.score_threshold,
                         png_dir=args.summary_png_dir, num_reads=args.num_reads)
        return

    if args.input is None or args.output_dir is None:
        parser.error("--input and --output-dir are required unless --combined-pdf is set")

    df = read_table(args.input)
    res_df = read_table(args.resources) if args.resources else None
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # secondary-accuracy metrics are plotted only when present in the table
    metrics = with_secondary_if_present(METRICS, df)

    if args.k_plot:
        # metric vs k, at a single score threshold; k-less tools as dashed lines
        score_filter = score_eq(args.score_threshold)
        for col_name, file_name in metrics:
            plot_xy(df, x="k", y=col_name,
                    output_path=args.output_dir / f"k_{file_name}.png",
                    filter=score_filter, hue="Base Name", baseline_hue="Tool",
                    baseline=True, title=col_name)
        # resource metrics (runtime, memory) as line vs k, like the accuracy metrics
        if res_df is not None:
            for y_col, file_name, title in RESOURCE_METRICS:
                plot_xy(res_df, x="k", y=y_col,
                        output_path=args.output_dir / f"{file_name}.png",
                        hue="Base Name", baseline_hue="Tool", baseline=True,
                        title=title)
    else:
        # metric vs score threshold, at a single k (resource metrics are k-only)
        k_filter = k_eq(args.k)
        for col_name, file_name in metrics:
            plot_xy(df, x="Score Threshold", y=col_name,
                    output_path=args.output_dir / f"{file_name}.png",
                    filter=k_filter, hue="Tool", baseline=False, title=col_name)


if __name__ == "__main__":
    main()
