#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run_snakemake(snakefile, config, cores=None, slurm=False, secondary_accuracy=False, dry_run=False):
    cmd = ["snakemake", "-s", snakefile, "--use-singularity", "--configfile", config]

    if slurm:
        cmd.extend(["--profile", "simx/profiles/slurm"])
    else:
        cmd.extend(["--cores", str(cores)])

    if secondary_accuracy:
        cmd.extend(["--config", "secondary_accuracy=True"])

    if dry_run:
        cmd.append("--dry-run")

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(result.returncode)


def main():
    parser = argparse.ArgumentParser(description="Run the aDNA benchmark pipeline")
    parser.add_argument("--test", action="store_true", help="Use test config (fruitfly, 5000 reads)")
    parser.add_argument("--simulate-only", action="store_true", help="Run simulation phase only")
    parser.add_argument("--benchmark-only", action="store_true", help="Run benchmark phase only")
    parser.add_argument("--cores", type=int, default=4, help="Number of cores (default: 4, ignored with --slurm)")
    parser.add_argument("--slurm", action="store_true", help="Use SLURM profile instead of local execution")
    parser.add_argument("--secondary-accuracy", action="store_true", help="Enable the (costly) secondary-alignment accuracy metric")
    parser.add_argument("-n", "--dry-run", action="store_true", help="Only show which jobs would run (snakemake --dry-run)")
    args = parser.parse_args()

    # Workflow paths are relative to the repo root.
    os.chdir(REPO_ROOT)
    config = "simx/config/config_test.yaml" if args.test else "simx/config/config.yaml"
    run_simulate = not args.benchmark_only
    run_benchmark = not args.simulate_only

    print(f"Using config: {config}")
    if args.slurm:
        print("Using SLURM profile")
    else:
        print(f"Cores: {args.cores}")

    if run_simulate:
        print("Running simulation...")
        run_snakemake("simx/Simulation.snake", config, args.cores, args.slurm, dry_run=args.dry_run)

    if run_benchmark:
        print("Running benchmark...")
        run_snakemake("simx/Benchmark.snake", config, args.cores, args.slurm,
                      secondary_accuracy=args.secondary_accuracy, dry_run=args.dry_run)

    print("Done!")


if __name__ == "__main__":
    main()
