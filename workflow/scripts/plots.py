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


def plot_metric(df, metric_col, output_path):
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df, x="Score Threshold", y=metric_col, hue="Tool", marker="o")
    plt.xlabel("Score Threshold")
    plt.ylabel(metric_col)
    plt.title(metric_col)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for col_name, file_name in METRICS:
        plot_metric(df, col_name, args.output_dir / f"{file_name}.png")


if __name__ == "__main__":
    main()
