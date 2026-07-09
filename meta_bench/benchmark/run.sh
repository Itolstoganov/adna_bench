#!/usr/bin/env bash
# Run the selected taxonomic profiler on the simulated samples.
#
#   ./run.sh                 run the profiler (config: meta_bench/benchmark/config.yaml)
#   ./run.sh -n              dry run (extra args pass through to snakemake)
#   CORES=40 ./run.sh
#   PROFILE=meta_bench/benchmark/profiles/slurm ./run.sh    (submit to SLURM)
#
# The profiler is chosen by `profiler:` in the config (aMeta-malt implemented).
# aMeta runs under Snakemake 7 with its own conda envs and SLURM profile; this
# wrapper only generates the samplesheet and dispatches. Score afterwards with
# ./score.sh.
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"
CONFIG="${CONFIG:-meta_bench/benchmark/config.yaml}"
CORES="${CORES:-20}"
PROFILE="${PROFILE:-}"

eval "$(python - "$CONFIG" <<'PY'
import sys, yaml
c = yaml.safe_load(open(sys.argv[1]))
print(f'PROFILER={c["profiler"]}')
print(f'AMETA_DIR={c["ameta_dir"]}')
print(f'AMETA_CONFIG={c["ameta_config"]}')
print(f'SIM_RESULTS={c["sim_results"]}')
print(f'RESULTS_DIR={c["results_dir"]}')
print('SAMPLES="{}"'.format(" ".join(c["samples"])))
PY
)"

RUNDIR="$RESULTS_DIR/$PROFILER"
mkdir -p "$RUNDIR/config"

# Samplesheet: point aMeta at the simulated reads (absolute paths).
{
  printf 'sample\tfastq\n'
  for s in $SAMPLES; do
    reads="$SIM_RESULTS/$s/reads.fastq.gz"
    [[ -f "$reads" ]] || { echo "missing reads: $reads (run the simulation first)" >&2; exit 1; }
    printf '%s\t%s\n' "$s" "$(readlink -f "$reads")"
  done
} > "$RUNDIR/config/samples.tsv"

case "$PROFILER" in
  aMeta-malt)
    cp "$AMETA_CONFIG" "$RUNDIR/config/config.yaml"
    profile_args=()
    [[ -n "$PROFILE" ]] && profile_args=(--profile "$(readlink -f "$PROFILE")")
    ( cd "$RUNDIR" && snakemake -s "$(readlink -f "$AMETA_DIR")/workflow/Snakefile" \
        --configfile config/config.yaml --use-conda --conda-frontend conda \
        "${profile_args[@]}" -j "$CORES" "$@" )
    ;;
  aMeta-ngslca)
    echo "profiler 'aMeta-ngslca' is not implemented yet (aligner + ngsLCA)." >&2
    exit 1 ;;
  *)
    echo "unknown profiler: '$PROFILER' (config: $CONFIG)" >&2; exit 1 ;;
esac

echo "Profiler '$PROFILER' finished. Score with: ./meta_bench/benchmark/score.sh"
