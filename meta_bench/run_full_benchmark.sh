#!/usr/bin/env bash
# Full metagenomic benchmark, end to end. Four steps:
#
#   1. databases   download the KrakenUniq and MALT databases
#   2. simulate    simulate 10 ancient metagenomic samples with gargammel
#   3. ameta       run aMeta (benchmark fork: KrakenUniq -> MALT + strobealign/ngsLCA
#                  -> authentication) with config/config_ngslca.yaml
#   4. summary     score the authentication results against the simulation ground truth,
#                  and summarise runtime and memory from the Snakemake benchmarks
#
# Usage (from anywhere; the script locates the repo itself):
#   meta_bench/run_full_benchmark.sh               # all four steps, in order
#   meta_bench/run_full_benchmark.sh 3 4           # only aMeta + summary
#   meta_bench/run_full_benchmark.sh simulate      # by name: databases|simulate|ameta|summary
#   meta_bench/run_full_benchmark.sh 2 3 -- -n     # args after -- go to snakemake (dry run)
#
# Env overrides:
#   CORES=40         cores for both snakemake runs (default 20)
#   RAM_GB=156       advisory memory budget for the aMeta scheduler (default 156)
#   N_SAMPLES=10     number of simulated samples (default 10)
#   CONDA=mamba      solver used to create missing envs (default conda)
# Step 1 needs ~800 GB free for KrakenUniq
# and ~390 GB for the MALT nt reference; the download targets are read out of the
# aMeta config, so they always match the paths aMeta will look in.
#
# This deliberately does not use run_ameta.sh (which adds a systemd cgroup hard RAM
# cap and picks a samplesheet); step 3 is the plain snakemake command.
set -euo pipefail

MB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"      # meta_bench/
REPO_ROOT="$(cd "$MB_DIR/.." && pwd)"
# Both snakemake workflows resolve their paths relative to the repo root.
cd "$REPO_ROOT"

RUNDIR="$MB_DIR/runs/sim_dataset"
AMETA_DIR="$MB_DIR/aMeta"
AMETA_CONFIG="$RUNDIR/config/config_ngslca.yaml"
SIM_RESULTS="$MB_DIR/results/sim"

CORES="${CORES:-20}"
RAM_GB="${RAM_GB:-156}"
N_SAMPLES="${N_SAMPLES:-10}"
CONDA="${CONDA:-conda}"
SIM_ENV="${SIM_ENV:-adna-eval}"
SIM_ENV_YML="$REPO_ROOT/workflow/envs/adna-eval.yml"
AMETA_ENV="${AMETA_ENV:-aMeta}"
AMETA_ENV_YML="$AMETA_DIR/workflow/envs/environment.yaml"

KU_URL_BASE="https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2022-06-16-STANDARD"
KU_TGZ="kuniq_standard_minus_kdb.20220616.tgz"
KU_KDB_SIZE=404451462724   # database.kdb is served separately from the tarball

# figshare file ids of the MALT reference + taxonomy maps (aMeta vignette,
# meta_bench/aMeta/vignettes/aMeta.qmd)
FIGSHARE_LIBRARY_FNA=37395181
FIGSHARE_SEQID2TAXID=38201982
FIGSHARE_ACC2TAXID=38201937

banner() { printf '\n=== step %s: %s ===\n' "$1" "$2" >&2; }

# Run a command inside a conda env (the two snakemake versions must not mix).
in_env() {
  local env="$1"; shift
  ( eval "$(conda shell.bash hook)"; conda activate "$env"; exec "$@" )
}

env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }

# Create a conda env from its yaml if it does not exist yet.
ensure_env() {
  local env="$1" yml="$2"
  if env_exists "$env"; then
    echo ">> conda env '$env' present" >&2
    return
  fi
  [[ -f "$yml" ]] || {
    echo "error: cannot create env '$env': $yml is missing" \
         "(aMeta is a submodule - run 'git submodule update --init meta_bench/aMeta')" >&2
    exit 1
  }
  echo ">> creating conda env '$env' from $yml (this takes a while)" >&2
  "$CONDA" env create -n "$env" -f "$yml"
}

# Build the gargammel fork which includes deaminated-position output. 
ensure_gargammel() {
  local src="$REPO_ROOT/gargammel"
  if [[ -x "$src/src/fragSim" && -x "$src/src/deamSim" ]]; then
    echo ">> gargammel fork already built" >&2
    return
  fi
  [[ -f "$src/Makefile" ]] || git -C "$REPO_ROOT" submodule update --init gargammel
  echo ">> building the gargammel fork in $src" >&2
  in_env "$SIM_ENV" bash -c "cd '$src' && make"
  [[ -x "$src/src/fragSim" && -x "$src/src/deamSim" ]] || {
    echo "error: the gargammel build did not produce src/fragSim + src/deamSim" >&2
    exit 1
  }
}

cfg_get() {
  local v
  v="$(sed -nE "s/^$1:[[:space:]]*//p" "$AMETA_CONFIG" | head -1)"
  v="${v%\"}"; v="${v#\"}"
  [[ -n "$v" ]] || { echo "error: key '$1' not found in $AMETA_CONFIG" >&2; exit 1; }
  printf '%s\n' "$v"
}

fetch() {
  local url="$1" out="$2"
  mkdir -p "$(dirname "$out")"
  echo ">> downloading $url -> $out" >&2
  curl -SL --fail --continue-at - \
       --retry 50 --retry-delay 15 --retry-all-errors \
       -o "$out" "$url"
}

fetch_figshare() { fetch "https://ndownloader.figshare.com/files/$1" "$2"; }

check_space() {
  local dir="$1" need_gb="$2" free_gb
  mkdir -p "$dir"
  free_gb=$(df -BG --output=avail "$dir" | tail -1 | tr -dc '0-9')
  if (( free_gb < need_gb )); then
    echo "WARNING: only ${free_gb}G free on $dir, ~${need_gb}G needed" >&2
  fi
}

##############################################################################
# 1. download the KrakenUniq and MALT databases
##############################################################################
step_databases() {
  banner 1 "download the KrakenUniq and MALT databases"

  local ku_db malt_fna seqid2taxid acc2taxid
  ku_db="$(cfg_get krakenuniq_db)"
  malt_fna="$(cfg_get malt_nt_fasta)"
  seqid2taxid="$(cfg_get malt_seqid2taxid_db)"
  acc2taxid="$(cfg_get malt_accession2taxid)"

  if [[ -s "$ku_db/taxDB" && -s "$ku_db/database.idx" ]]; then
    echo ">> KrakenUniq index present, skipping" >&2
  else
    check_space "$ku_db" 20
    fetch "$KU_URL_BASE/$KU_TGZ" "$ku_db/$KU_TGZ"
    tar -xzf "$ku_db/$KU_TGZ" -C "$ku_db"
  fi
  if [[ "$(stat -c %s "$ku_db/database.kdb" 2>/dev/null || echo 0)" == "$KU_KDB_SIZE" ]]; then
    echo ">> KrakenUniq database.kdb present, skipping" >&2
  else
    check_space "$ku_db" 380
    fetch "$KU_URL_BASE/database.kdb" "$ku_db/database.kdb"
  fi

  # MALT: nt reference FASTA + the two taxonomy maps aMeta builds its DB from
  # (library.fna.fai is created later by the workflow itself).
  if [[ -s "$malt_fna" ]]; then
    echo ">> $(basename "$malt_fna") present, skipping" >&2
  else
    check_space "$(dirname "$malt_fna")" 500   # 105 GB gzipped + 348 GB unpacked
    fetch_figshare "$FIGSHARE_LIBRARY_FNA" "$malt_fna.gz"
    gunzip "$malt_fna.gz"
  fi
  [[ -s "$seqid2taxid" ]] || fetch_figshare "$FIGSHARE_SEQID2TAXID" "$seqid2taxid"
  [[ -s "$acc2taxid" ]]   || fetch_figshare "$FIGSHARE_ACC2TAXID" "$acc2taxid"
}

##############################################################################
# 2. simulate the samples
##############################################################################
step_simulate() {
  banner 2 "simulate $N_SAMPLES ancient metagenomic samples with gargammel"

  ensure_gargammel

  in_env "$SIM_ENV" snakemake -s "$MB_DIR/workflow/Simulation.snake" \
    --cores "$CORES" --config "n_samples=$N_SAMPLES" --rerun-incomplete "${sm_args[@]}"

  # aMeta reads the reads through the run dir's fastq/ symlinks (config/*.tsv).
  mkdir -p "$RUNDIR/fastq"
  for i in $(seq 1 "$N_SAMPLES"); do
    ln -sfn "$SIM_RESULTS/sample$i/reads.fastq.gz" "$RUNDIR/fastq/sample$i.fastq.gz"
  done
}

##############################################################################
# 3. run aMeta
##############################################################################
step_ameta() {
  banner 3 "run aMeta (KrakenUniq -> MALT + strobealign/ngsLCA -> authentication)"

  local total_mb=$((RAM_GB * 1024)) default_mb=4000
  # Advisory memory costs of the heavy rules, so the scheduler does not start two
  # of them at once; clamped so a small RAM_GB still leaves every job schedulable.
  clamp() { local v=$1 cap=$((total_mb - default_mb)); (( v > cap )) && v=$cap; echo "$v"; }
  local ku_mb malt_build_mb malt_mb malt_extract_mb
  ku_mb="$(clamp 96000)"            # krakenuniq_preload_size 80G + overhead
  malt_build_mb="$(clamp 147456)"   # malt_max_heap 128G + overhead
  malt_mb="$(clamp 147456)"
  malt_extract_mb="$(clamp 40000)"  # MaltExtract -Xmx32G + overhead

  # Run from the run dir: aMeta's common.smk also loads a hardcoded
  # `configfile: "config/config.yaml"` relative to cwd and deep-merges it under
  # --configfile, so the run dir must be the one holding both configs.
  ( cd "$RUNDIR" && in_env "$AMETA_ENV" snakemake \
      --snakefile "$AMETA_DIR/workflow/Snakefile" \
      --configfile config/config_ngslca.yaml \
      --use-conda --conda-frontend conda \
      --cores "$CORES" \
      --resources "mem_mb=$total_mb" \
      --default-resources "mem_mb=$default_mb" \
      --set-resources \
        "KrakenUniq:mem_mb=$ku_mb" \
        "Build_Malt_DB:mem_mb=$malt_build_mb" \
        "Malt:mem_mb=$malt_mb" \
        "Malt_Extract:mem_mb=$malt_extract_mb" \
      --rerun-incomplete --printshellcmds "${sm_args[@]}" )
}

##############################################################################
# 4. summarise against the ground truth
##############################################################################
step_summary() {
  banner 4 "score the authentication results against the simulation ground truth"

  in_env "$SIM_ENV" python "$MB_DIR/scripts/auth_summary.py" \
    --results "$RUNDIR/results" \
    --truth "$SIM_RESULTS" \
    --outdir "$RUNDIR/results/meta_summary"

  in_env "$SIM_ENV" python "$MB_DIR/scripts/bench_summary.py" \
    --rundir "$RUNDIR" \
    --outdir "$RUNDIR/results/meta_summary/bench"
}

##############################################################################
# main
##############################################################################
steps=()
args=()
sm_args=()          # everything after `--` is passed on to the snakemake runs
while [[ $# -gt 0 ]]; do
  if [[ "$1" == "--" ]]; then shift; sm_args=("$@"); break; fi
  args+=("$1"); shift
done

if [[ ${#args[@]} -eq 0 ]]; then
  steps=(1 2 3 4)
else
  for a in "${args[@]}"; do
    case "$a" in
      1|databases) steps+=(1) ;;
      2|simulate)  steps+=(2) ;;
      3|ameta)     steps+=(3) ;;
      4|summary)   steps+=(4) ;;
      -h|--help)   awk 'NR>1{ if (/^#/) print; else exit }' "${BASH_SOURCE[0]}"; exit 0 ;;
      *) echo "usage: $0 [1|databases] [2|simulate] [3|ameta] [4|summary]" >&2; exit 1 ;;
    esac
  done
fi

# Conda envs: create whichever the selected steps need and do not exist yet.
steps_sorted="$(printf '%s\n' "${steps[@]}" | sort -u)"
grep -qE '^(2|4)$' <<<"$steps_sorted" && ensure_env "$SIM_ENV" "$SIM_ENV_YML"
grep -qE '^3$'     <<<"$steps_sorted" && ensure_env "$AMETA_ENV" "$AMETA_ENV_YML"

for s in $steps_sorted; do
  case "$s" in
    1) step_databases ;;
    2) step_simulate ;;
    3) step_ameta ;;
    4) step_summary ;;
  esac
done

echo
echo "Done. Results: $RUNDIR/results (meta_summary/ holds the ground-truth heatmap and metrics, meta_summary/bench/ the runtime/memory summary)."
