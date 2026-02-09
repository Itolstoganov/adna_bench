#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

METRICS = [
    ("% Mapped endogenous reads", "mapped_endogenous"),
    ("% Accuracy (endogenous only)", "accuracy"),
    ("% Mapped bacterial reads", "mapped_bacterial"),
    ("% Mapped contaminated reads", "mapped_contaminated"),
]


def plot_metric(df, metric_col, output_path, k):
    filtered_df = df[(df["k"] == k) | (df["k"].astype(str).str.strip() == "nan")]
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=filtered_df, x="Score Threshold", y=metric_col, hue="Tool", marker="o")
    plt.xlabel("Score Threshold")
    plt.ylabel(metric_col)
    plt.title(metric_col)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def plot_k_metric(df, metric_col, output_path, score_threshold):
    # print(df["Base Name"])
    # print(df["k"].astype(str).str.strip())
    k_df = df[(df["k"].astype(str).str.strip() != "nan")].copy()
    # print(k_df)
    k_df["k"] = k_df["k"].astype(int)
    k_df = k_df[k_df["Score Threshold"] == score_threshold]
    if k_df.empty:
        return

    plt.figure(figsize=(10, 6))
    sns.lineplot(data=k_df, x="k", y=metric_col, hue="Base Name", marker="o")
    plt.xlabel("k")
    plt.ylabel(metric_col)
    plt.title(f"{metric_col} (score threshold = {score_threshold})")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(output_path)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--k-plot", action="store_true")
    parser.add_argument("--score-threshold", type=int, default=60)
    parser.add_argument("--k", type=int)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.k_plot:
        for col_name, file_name in METRICS:
            plot_k_metric(df, col_name, args.output_dir / f"k_{file_name}.pdf", args.score_threshold)
    else:
        for col_name, file_name in METRICS:
            plot_metric(df, col_name, args.output_dir / f"{file_name}.pdf", args.k)


if __name__ == "__main__":
    main()
