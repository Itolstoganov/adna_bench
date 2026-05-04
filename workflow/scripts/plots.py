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

COMBINED_K_METRICS = [
    ("% Accuracy (endogenous only)", "accuracy"),
    ("% Accuracy (deaminated endogenous only)", "accuracy_deam"),
    ("% Mapped endogenous reads", "mapped_endogenous"),
    ("% Mapped deaminated endogenous reads", "mapped_endogenous_deam"),
]


def plot_metric(df, metric_col, output_path, k):
    filtered_df = df[(df["k"] == k) | (df["k"].astype(str).str.strip() == "nan")]
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=filtered_df, x="Score Threshold", y=metric_col, hue="Tool", marker="o")
    plt.xticks(sorted(filtered_df["Score Threshold"].unique()))
    plt.xlabel("Score Threshold")
    plt.ylabel(metric_col)
    plt.title(metric_col)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def render_k_metric(df, metric_col, score_threshold, title=None):
    k_df = df[(df["k"].astype(str).str.strip() != "nan")].copy()
    if not k_df.empty:
        k_df["k"] = k_df["k"].astype(int)
        k_df = k_df[k_df["Score Threshold"] == score_threshold]

    non_k_df = df[(df["k"].astype(str).str.strip() == "nan")].copy()
    non_k_df = non_k_df[non_k_df["Score Threshold"] == score_threshold]

    if k_df.empty and non_k_df.empty:
        return False

    all_names = list(k_df["Base Name"].unique()) + list(non_k_df["Tool"].unique())
    palette = dict(zip(all_names, sns.color_palette("tab10", n_colors=len(all_names))))

    plt.figure(figsize=(10, 6))
    if not k_df.empty:
        sns.lineplot(data=k_df, x="k", y=metric_col, hue="Base Name", palette=palette, marker="o")

    for _, row in non_k_df.iterrows():
        plt.axhline(y=row[metric_col], linestyle="--", color=palette[row["Tool"]], label=row["Tool"])

    if not k_df.empty:
        plt.xticks(sorted(k_df["k"].unique()))
    plt.xlabel("k")
    plt.ylabel(metric_col)
    plt.title(title if title is not None else metric_col)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    return True


def plot_k_metric(df, metric_col, output_path, score_threshold):
    if render_k_metric(df, metric_col, score_threshold):
        plt.savefig(output_path, dpi=150)
    plt.close()


def _build_combined_frame(dataset_inputs, score_threshold):
    frames = []
    for _, path, row_label, col_label in dataset_inputs:
        df = pd.read_csv(path, sep="\t")
        df = df[df["Score Threshold"] == score_threshold].copy()
        df["dataset_row"] = row_label
        df["dataset_col"] = col_label
        frames.append(df)
    if not frames:
        return None
    df = pd.concat(frames, ignore_index=True)
    is_k_aware = df["k"].astype(str).str.strip() != "nan"
    df["k_aware"] = is_k_aware
    df["Display"] = df["Base Name"].where(is_k_aware, df["Tool"])
    return df


def plot_combined_k_pdf(dataset_inputs, output_path, score_threshold):
    df = _build_combined_frame(dataset_inputs, score_threshold)
    if df is None or df.empty:
        return

    k_aware_df = df[df["k_aware"]].copy()
    if k_aware_df.empty:
        return
    k_aware_df["k"] = k_aware_df["k"].astype(int)
    k_values = sorted(k_aware_df["k"].unique())

    non_k_df = df[~df["k_aware"]].copy()
    if not non_k_df.empty and k_values:
        broadcast = non_k_df.loc[non_k_df.index.repeat(len(k_values))].copy()
        broadcast["k"] = k_values * len(non_k_df)
        plot_df = pd.concat([k_aware_df, broadcast], ignore_index=True)
    else:
        plot_df = k_aware_df

    seen = set()
    display_order = []
    for name in list(k_aware_df["Display"].unique()) + list(non_k_df["Display"].unique()):
        if name not in seen:
            seen.add(name)
            display_order.append(name)
    palette = dict(zip(display_order, sns.color_palette("tab10", n_colors=len(display_order))))
    non_k_displays = set(non_k_df["Display"].unique())

    row_order = sorted(plot_df["dataset_row"].unique())
    col_order = sorted(plot_df["dataset_col"].unique())

    with PdfPages(output_path) as pdf:
        for metric_col, _ in COMBINED_K_METRICS:
            g = sns.relplot(
                data=plot_df,
                x="k", y=metric_col,
                hue="Display", hue_order=display_order, palette=palette,
                style="k_aware", style_order=[True, False],
                row="dataset_row", col="dataset_col",
                row_order=row_order, col_order=col_order,
                kind="line", marker="o",
                facet_kws={"sharey": False, "legend_out": True},
                legend=False,
            )
            g.set_axis_labels("k", metric_col)
            g.set_titles(row_template="{row_name}", col_template="{col_name}")
            for ax in g.axes.flat:
                ax.set_xticks(k_values)

            handles = []
            for display in display_order:
                is_k_aware = display not in non_k_displays
                handles.append(Line2D(
                    [], [],
                    color=palette[display],
                    linestyle="-" if is_k_aware else "--",
                    marker="o" if is_k_aware else "",
                ))
            g.figure.legend(
                handles, display_order, title="Tool",
                loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False,
            )

            g.figure.suptitle(metric_col)
            g.tight_layout()
            pdf.savefig(g.figure, bbox_inches="tight")
            plt.close(g.figure)


def plot_runtime(res_df, output_path, k):
    filtered = res_df[(res_df["k"] == k) | (res_df["k"].astype(str).str.strip() == "nan")]
    if filtered.empty:
        return

    plt.figure(figsize=(10, 6))
    sns.barplot(data=filtered, x="Tool", y="Mapping Time (s)")
    plt.xlabel("Tool")
    plt.ylabel("Mapping Time (s)")
    plt.title("Mapping Time")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def plot_user_time(res_df, output_path, k):
    filtered = res_df[(res_df["k"] == k) | (res_df["k"].astype(str).str.strip() == "nan")]
    if filtered.empty:
        return

    plt.figure(figsize=(10, 6))
    sns.barplot(data=filtered, x="Tool", y="User Time (s)")
    plt.xlabel("Tool")
    plt.ylabel("User Time (s)")
    plt.title("User Time")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def plot_k_user_time(res_df, output_path):
    k_df = res_df[res_df["k"].astype(str).str.strip() != "nan"].copy()
    k_df["k"] = k_df["k"].astype(int)

    non_k_df = res_df[res_df["k"].astype(str).str.strip() == "nan"].copy()

    if k_df.empty and non_k_df.empty:
        return

    all_names = list(k_df["Base Name"].unique()) + list(non_k_df["Tool"].unique())
    palette = dict(zip(all_names, sns.color_palette("tab10", n_colors=len(all_names))))

    plt.figure(figsize=(10, 6))
    if not k_df.empty:
        sns.lineplot(data=k_df, x="k", y="User Time (s)", hue="Base Name", palette=palette, marker="o")

    for _, row in non_k_df.iterrows():
        plt.axhline(y=row["User Time (s)"], linestyle="--", color=palette[row["Tool"]], label=row["Tool"])

    if not k_df.empty:
        plt.xticks(sorted(k_df["k"].unique()))
    plt.xlabel("k")
    plt.ylabel("User Time (s)")
    plt.title("User Time")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def plot_k_runtime(res_df, output_path):
    k_df = res_df[res_df["k"].astype(str).str.strip() != "nan"].copy()
    k_df["k"] = k_df["k"].astype(int)

    non_k_df = res_df[res_df["k"].astype(str).str.strip() == "nan"].copy()

    if k_df.empty and non_k_df.empty:
        return

    all_names = list(k_df["Base Name"].unique()) + list(non_k_df["Tool"].unique())
    palette = dict(zip(all_names, sns.color_palette("tab10", n_colors=len(all_names))))

    plt.figure(figsize=(10, 6))
    if not k_df.empty:
        sns.lineplot(data=k_df, x="k", y="Mapping Time (s)", hue="Base Name", palette=palette, marker="o")

    for _, row in non_k_df.iterrows():
        plt.axhline(y=row["Mapping Time (s)"], linestyle="--", color=palette[row["Tool"]], label=row["Tool"])

    if not k_df.empty:
        plt.xticks(sorted(k_df["k"].unique()))
    plt.xlabel("k")
    plt.ylabel("Mapping Time (s)")
    plt.title("Mapping Time")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


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
    parser.add_argument("--dataset", nargs=4, action="append", default=[],
                        metavar=("NAME", "PATH", "ROW", "COL"),
                        help="Dataset name, accuracy table path, row label, "
                             "column label (repeatable, used with --combined-pdf).")
    args = parser.parse_args()

    if args.combined_pdf:
        if not args.dataset:
            parser.error("--combined-pdf requires at least one --dataset entry")
        args.combined_pdf.parent.mkdir(parents=True, exist_ok=True)
        plot_combined_k_pdf(args.dataset, args.combined_pdf, args.score_threshold)
        return

    if args.input is None or args.output_dir is None:
        parser.error("--input and --output-dir are required unless --combined-pdf is set")

    df = pd.read_csv(args.input, sep="\t")
    res_df = pd.read_csv(args.resources, sep="\t") if args.resources else None
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.k_plot:
        for col_name, file_name in METRICS:
            plot_k_metric(df, col_name, args.output_dir / f"k_{file_name}.png", args.score_threshold)
        if res_df is not None:
            plot_k_runtime(res_df, args.output_dir / "k_runtime.png")
            plot_k_user_time(res_df, args.output_dir / "k_user_time.png")
    else:
        for col_name, file_name in METRICS:
            plot_metric(df, col_name, args.output_dir / f"{file_name}.png", args.k)
        if res_df is not None:
            plot_runtime(res_df, args.output_dir / "runtime.png", args.k)
            plot_user_time(res_df, args.output_dir / "user_time.png", args.k)


if __name__ == "__main__":
    main()
