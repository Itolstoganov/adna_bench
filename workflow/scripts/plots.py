#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from pathlib import Path

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

def _palette(names):
    """Stable tab10 color mapping for a list of (possibly duplicated) names."""
    seen = []
    for n in names:
        if n not in seen:
            seen.append(n)
    return dict(zip(seen, sns.color_palette("tab10", n_colors=len(seen))))


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

    plt.figure(figsize=(10, 6))
    if not varying.empty:
        sns.lineplot(data=varying, x=x, y=y,
                     hue=hue if hue else None, palette=palette, marker="o")
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
    base_displays = set(base_df[hue].unique())

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
        style="k_aware", style_order=[True, False],
        kind="line", marker="o",
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
            marker="o" if k_aware else "",
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
        df = pd.read_csv(path, sep="\t")
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
        df = pd.read_csv(path, sep="\t")
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
SUMMARY_DELTA_BASELINE = "bwaaln"

# Friendlier metric labels for summary plot titles / axes. Keys are the table
# column names; everything not listed is shown verbatim.
SUMMARY_METRIC_LABELS = {
    "% Accuracy (endogenous only)": "% Accuracy",
    "% Accuracy (deaminated endogenous only)": "% Accuracy (deaminated reads only)",
}


def _summary_label(col):
    return SUMMARY_METRIC_LABELS.get(col, col)


def _apply_delta(df, y, baseline):
    """Subtract `baseline`'s per-dataset value of `y` from every row.

    Returns a copy with `y` replaced by the difference to the baseline tool. Rows
    of datasets where the baseline is absent get NaN and are dropped downstream.
    """
    base = (df[df["Display"] == baseline][["dataset", y]]
            .rename(columns={y: "__base"}))
    merged = df.merge(base, on="dataset", how="left").copy()
    merged[y] = merged[y] - merged["__base"]
    return merged.drop(columns="__base")


def plot_summary_metric(df, y, dataset_order, *, pdf=None, png_path=None,
                        title=None, delta_baseline=None):
    """Point/line plot of metric `y` with datasets on the x axis, one line per tool.

    When `delta_baseline` is set, values are plotted as the difference to that
    tool (a dashed line marks zero). The figure is written to `pdf` and/or
    `png_path` (both may be given to emit the same plot to each).
    """
    if y not in df.columns:
        return
    if delta_baseline is not None:
        df = _apply_delta(df, y, delta_baseline)
    sub = df[df[y].notna()]
    if sub.empty:
        return
    names = list(dict.fromkeys(sub["Display"]))
    palette = _palette(names)
    plt.figure(figsize=(max(10, len(dataset_order) * 1.2), 6))
    if delta_baseline is not None:
        plt.axhline(0, color="gray", linewidth=0.8, linestyle="--")
    sns.pointplot(data=sub, x="dataset", y=y, hue="Display",
                  order=dataset_order, hue_order=names, palette=palette,
                  errorbar=None)
    label = _summary_label(y)
    plt.xlabel("Dataset")
    plt.ylabel(f"Δ {label} (vs {delta_baseline})" if delta_baseline is not None else label)
    plt.title(title if title is not None else y)
    plt.xticks(rotation=45, ha="right")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    if pdf is not None:
        pdf.savefig(bbox_inches="tight")
    if png_path is not None:
        plt.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_summary_pdf(summary_inputs, output_path, k, score_threshold, png_dir=None):
    # summary_inputs: list of (NAME, ACC_PATH, RES_PATH, GENOME, READ_LEN, VARIATION)
    short_names = _build_short_names(
        [(genome, read_len, var) for _n, _a, _r, genome, read_len, var in summary_inputs]
    )
    dataset_order = list(dict.fromkeys(short_names))
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
            for y_col, slug, title in RESOURCE_METRICS:
                plot_summary_metric(res_df, y_col, dataset_order, pdf=pdf,
                                    png_path=png_path(slug), title=title)


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
    parser.add_argument("--summary-dataset", nargs=6, action="append", default=[],
                        metavar=("NAME", "ACC_PATH", "RES_PATH", "GENOME",
                                 "READ_LEN", "VARIATION"),
                        help="Dataset name, accuracy table path, resource table "
                             "path, genome, read length, variation (repeatable, "
                             "used with --summary-pdf).")
    args = parser.parse_args()

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
                         png_dir=args.summary_png_dir)
        return

    if args.input is None or args.output_dir is None:
        parser.error("--input and --output-dir are required unless --combined-pdf is set")

    df = pd.read_csv(args.input, sep="\t")
    res_df = pd.read_csv(args.resources, sep="\t") if args.resources else None
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
