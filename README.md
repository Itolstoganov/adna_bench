# adna_bench

Benchmarks for ancient-DNA tools on gargammel-simulated reads with per-read ground truth.

- `simx/` — read mapping accuracy benchmark (strobealign, bwa aln/mem, bowtie2, SAFARI).
- `meta_bench/` — metagenomic taxonomic profiler benchmark (aMeta).

All commands are run from the repository root.

## Install

Requirements: conda (or mamba), git, and for `simx` also `cargo` (to build strobealign), `curl`, `wget` and `unzip`.

```bash
git clone https://github.com/Itolstoganov/adna_bench.git adna_bench && cd adna_bench
./install.sh            # submodules, adna-eval conda env
conda activate adna-eval
```

## Aligner accuracy benchmark (simx)

1. Copy the example config and adjust datasets, aligner versions and parameters:

   ```bash
   cp simx/config/config_example.yaml simx/config/config.yaml
   ```

2. Run simulation and benchmark:

   ```bash
   simx/run.py --cores 32
   # SLURM (set account/partition in simx/profiles/slurm/config.yaml first):
   simx/run.py --slurm
   ```

   `--simulate-only` / `--benchmark-only` run one stage
   `-n` performs a dry run.

Results are written to `simx/result/`: per-dataset `accuracy_table.tsv`, `resource_table.tsv` and `plots/`, plus `summary_plots.pdf`.

## Metagenomic profiling benchmark (meta_bench)

1. Simulate the samples (composition and parameters in `meta_bench/config/config.yaml`; genomes are downloaded per `meta_bench/config/accessions.tsv`):

   ```bash
   snakemake -s meta_bench/workflow/Simulation.snake --cores 8
   # SLURM: add --profile meta_bench/benchmark/profiles/slurm -j 100
   ```

   Per-sample reads and per-read ground truth are written to `meta_bench/results/sim/{sample}/` (`reads.fastq.gz`, `truth.tsv`).

2. Set the KrakenUniq, MALT and NCBI taxonomy database paths in `meta_bench/benchmark/aMeta_config.yaml`, then run aMeta to obtain authentication scores:

   ```bash
   conda env create -n aMeta -f meta_bench/aMeta/workflow/envs/environment.yaml
   conda activate aMeta
   CORES=40 meta_bench/benchmark/run.sh     # SLURM: PROFILE=meta_bench/benchmark/profiles/slurm
   conda activate adna-eval
   meta_bench/benchmark/score.sh
   ```

   The profiler is selected by `profiler:` in `meta_bench/benchmark/config.yaml`. Scores are written to `meta_bench/results/benchmark/<profiler>/scores/`.
