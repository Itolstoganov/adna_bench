#!/usr/bin/env bash
# Adapt the profiler output to the normalized contract and score it against the
# simulation ground truth. Run after ./run.sh has produced the profiler results.
#
#   ./score.sh        (config: meta_bench/benchmark/config.yaml)
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"
CONFIG="${CONFIG:-meta_bench/benchmark/config.yaml}"

eval "$(python - "$CONFIG" <<'PY'
import sys, yaml
c = yaml.safe_load(open(sys.argv[1]))
print(f'PROFILER={c["profiler"]}')
print(f'SIM_RESULTS={c["sim_results"]}')
print(f'RESULTS_DIR={c["results_dir"]}')
print(f'DETECT_THRESHOLD={c.get("detect_threshold", 1)}')
print(f'MIN_TRUTH_READS={c.get("min_truth_reads", 1)}')
print('SAMPLES="{}"'.format(" ".join(c["samples"])))
PY
)"

RUNDIR="$RESULTS_DIR/$PROFILER"
BENCH_SCRIPTS="meta_bench/benchmark/scripts"
DETECTIONS="$RUNDIR/taxon_abundance.tsv"

# 1. Adapter: profiler-native output -> normalized sample/taxid/taxon/count.
case "$PROFILER" in
  aMeta-malt)
    python "$BENCH_SCRIPTS/adapter_malt.py" \
      --matrix "$RUNDIR/results/MALT_ABUNDANCE_MATRIX_SAM/malt_abundance_matrix_sam.txt" \
      --out "$DETECTIONS" ;;
  aMeta-ngslca)
    python "$BENCH_SCRIPTS/adapter_ngslca.py" ;;
  *)
    echo "unknown profiler: '$PROFILER'" >&2; exit 1 ;;
esac

# 2. Score against the per-read simulation truth.
python "$BENCH_SCRIPTS/score.py" \
  --sim-results "$SIM_RESULTS" --samples $SAMPLES \
  --detections "$DETECTIONS" --out "$RUNDIR/scores" \
  --detect-threshold "$DETECT_THRESHOLD" --min-truth-reads "$MIN_TRUTH_READS"
