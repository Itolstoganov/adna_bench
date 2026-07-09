# meta_bench — metagenomic aDNA simulation & benchmark

Simulate ancient-DNA metagenomic datasets (with per-read ground truth) and
benchmark taxonomic profilers against them. Two stages:

1. **Simulation** (`meta_bench/workflow/`, Snakemake 9) — reproduces the simulated dataset from [Pochon et al., 2023](https://doi.org/10.1186/s13059-023-03083-9), with additional  per-read ground truth table (taxon + coordinates + deaminated positions).
2. **Benchmark** (`meta_bench/benchmark/`) — runs a profiler on the simulated reads
   and scores it against the truth. 

## Requirements

- **conda/mamba**.
- **git** with submodule support.
- Build toolchain for gargammel: `make`, C/C++ compiler, `cmake`, `zlib`, `gsl`,
  `seq-gen`, `perl` — all provided by the `adna-eval` conda env
  (`workflow/envs/adna-eval.yml`). 

## Install

```bash
git clone <adna_bench-url> && cd adna_bench
./install.sh                 # submodules + adna-eval env + build gargammel fork
conda activate adna-eval
```

`install.sh` is idempotent. Use `CONDA=mamba ./install.sh` if you have mamba.

## Run the simulation

Reference genomes are pinned in `meta_bench/config/accessions.tsv` (RefSeq
accessions; Y. pestis = CO92). One row (`Vermamoeba_vermiformis`) is left blank and
resolved by name at download time; run
`meta_bench/workflow/scripts/resolve_accessions.py` to pin any remaining rows.

```bash
# from the repo root, env active
snakemake -s meta_bench/workflow/Simulation.snake --cores 8
# SLURM: snakemake -s meta_bench/workflow/Simulation.snake \
#          --profile meta_bench/benchmark/profiles/slurm -j 100
```

Outputs per sample under `meta_bench/results/sim/{sample}/`:
- `reads.fastq.gz` — single-end simulated reads (ancient + modern passes).
- `truth.tsv` — per read: `read_id class taxon taxid ref_id start end strand n_deam pass`.
- `{ancient,modern}/{ground_truth.bed, damage.bed, truth.tsv}` — per-pass detail.

Composition/params live in `meta_bench/config/config.yaml` (samples, seed, read
length, damage, per-sample bact/cont/endo proportions, human genome).

## Run the benchmark

Needs the aMeta submodule and its databases (KrakenUniq DB, MALT reference, NCBI
taxonomy). Set the DB paths in `meta_bench/benchmark/aMeta_config.yaml` (same paths
as `meta_bench/runs/sim_dataset/config/config.yaml`), and do aMeta's one-time env
setup (`meta_bench/runs/sim_dataset/setup_envs.sh`).

```bash
# 1. profile the simulated samples (aMeta under its own SM7 conda envs)
CORES=40 meta_bench/benchmark/run.sh
#    SLURM: PROFILE=meta_bench/benchmark/profiles/slurm meta_bench/benchmark/run.sh

# 2. adapt profiler output -> normalized table, then score vs. ground truth
meta_bench/benchmark/score.sh
```

Scores land in `meta_bench/results/benchmark/aMeta-malt/scores/`
(`per_sample_scores.tsv`, `summary.tsv`, `calls.tsv`). The profiler backend is set
by `profiler:` in `meta_bench/benchmark/config.yaml`. Bowtie2 + mapDamage are
disabled in the aMeta config (the benchmark scores the MALT profile, not the
pathogen screen).

## Notes

- Run all `snakemake` commands from the repo root (paths are repo-root-relative).
- The simulation (Snakemake 9) and aMeta (Snakemake 7) are intentionally separate
  environments; do not mix them.
- `meta_bench/scripts/aMeta/` (upstream reference scripts) is a separate clone and
  not required — the composition tables it provided are copied into
  `meta_bench/config/`.
