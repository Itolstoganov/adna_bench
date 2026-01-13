#!/usr/bin/env python3

import argparse
import subprocess
import sys


def run_snakemake(snakefile, config, cores=None, slurm=False):
    cmd = ["snakemake", "-s", snakefile, "--use-singularity", "--configfile", config]

    if slurm:
        cmd.extend(["--profile", "slurm"])
    else:
        cmd.extend(["--cores", str(cores)])

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(result.returncode)


def main():
    parser = argparse.ArgumentParser(description="Run the aDNA benchmark pipeline")
    parser.add_argument("--test", action="store_true", help="Use test config (ecoli, 5000 reads)")
    parser.add_argument("--simulate-only", action="store_true", help="Run simulation phase only")
    parser.add_argument("--benchmark-only", action="store_true", help="Run benchmark phase only")
    parser.add_argument("--cores", type=int, default=4, help="Number of cores (default: 4, ignored with --slurm)")
    parser.add_argument("--slurm", action="store_true", help="Use SLURM profile instead of local execution")
    args = parser.parse_args()

    config = "configs/config_test.yaml" if args.test else "configs/config.yaml"
    run_simulate = not args.benchmark_only
    run_benchmark = not args.simulate_only

    print(f"Using config: {config}")
    if args.slurm:
        print("Using SLURM profile")
    else:
        print(f"Cores: {args.cores}")

    if run_simulate:
        print("Running simulation...")
        run_snakemake("workflow/Simulation.snake", config, args.cores, args.slurm)

    if run_benchmark:
        print("Running benchmark...")
        run_snakemake("workflow/Benchmark.snake", config, args.cores, args.slurm)

    print("Done!")


if __name__ == "__main__":
    main()
