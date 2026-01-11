#!/bin/bash
set -e

CONFIG="configs/config.yaml"
CORES=8
RUN_SIMULATE=true
RUN_BENCHMARK=true

usage() {
    echo "Usage: $0 [--test] [--simulate-only] [--benchmark-only] [--cores N]"
    echo ""
    echo "Options:"
    echo "  --test            Use test config (ecoli, 5000 reads)"
    echo "  --simulate-only   Run simulation phase only"
    echo "  --benchmark-only  Run benchmark phase only"
    echo "  --cores N         Number of cores (default: 8)"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --test)
            CONFIG="configs/config_test.yaml"
            shift
            ;;
        --simulate-only)
            RUN_BENCHMARK=false
            shift
            ;;
        --benchmark-only)
            RUN_SIMULATE=false
            shift
            ;;
        --cores)
            CORES="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

echo "Using config: $CONFIG"
echo "Cores: $CORES"

if $RUN_SIMULATE; then
    echo "Running simulation..."
    snakemake -s workflow/Simulation.snake --use-singularity --configfile "$CONFIG" --cores "$CORES"
fi

if $RUN_BENCHMARK; then
    echo "Running benchmark..."
    snakemake -s workflow/Benchmark.snake --use-singularity --configfile "$CONFIG" --cores "$CORES"
fi

echo "Done!"
