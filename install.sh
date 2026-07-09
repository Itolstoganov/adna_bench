#!/usr/bin/env bash
#   git clone <adna_bench> && cd adna_bench && ./install.sh
#
# Env overrides:  ENV_NAME (default adna-eval)   CONDA (default conda; use mamba if you have it)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_ROOT"

ENV_NAME="${ENV_NAME:-adna-eval}"
CONDA="${CONDA:-conda}"

echo "[1/3] Submodules ------------------------------------------------------"
# gargammel fork which outputs deaminated positions. Do not use the conda `gargammel`, which lacks that feature.
git submodule update --init gargammel
# aMeta (optional; only for the metagenomic benchmark).
git submodule update --init meta_bench/aMeta \
  || echo "  WARN: aMeta submodule init failed (needs GitHub SSH access)."

echo "[2/3] Conda env '$ENV_NAME' -----------------------------------------"
if $CONDA env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "  exists; updating from workflow/envs/adna-eval.yml"
  $CONDA env update -n "$ENV_NAME" -f workflow/envs/adna-eval.yml
else
  $CONDA env create -n "$ENV_NAME" -f workflow/envs/adna-eval.yml
fi

echo "[3/3] Building the gargammel fork (make) ----------------------------"
$CONDA run -n "$ENV_NAME" bash -c "cd gargammel && make"
if [[ -x gargammel/src/fragSim && -x gargammel/src/deamSim ]]; then
  echo "  gargammel built OK (src/fragSim, src/deamSim present)."
else
  echo "  ERROR: gargammel build did not produce src/fragSim + src/deamSim." >&2
  exit 1
fi

cat <<EOF

Done. Next:
  conda activate $ENV_NAME

  # Simulate the metagenomic datasets (writes meta_bench/results/sim/):
  snakemake -s meta_bench/workflow/Simulation.snake --cores 8

  # Benchmark (needs aMeta + its databases): see meta_bench/README.md
EOF
