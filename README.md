# adna_bench

Benchmarks for ancient-DNA tools on [gargammel](https://github.com/grenaud/gargammel)-simulated reads with per-read ground truth.

- `simx/` — read mapping accuracy benchmark ([strobealign](https://github.com/ksahlin/strobealign), [bwa](https://github.com/lh3/bwa) aln/mem, [bowtie2](https://github.com/BenLangmead/bowtie2)).
- `meta_bench/` — metagenomic taxonomic profiler benchmark ([aMeta](https://github.com/NBISweden/aMeta)).

All commands are run from the repository root.

## Install

Requirements: conda (or mamba), and for `simx` also `cargo` (to build strobealign), `curl`, `wget` and `unzip`.

`simx` additionally requires Singularity/Apptainer, with permission to pull and run containers from Docker Hub
(the g2gtools steps of the simulation run in `docker://churchilllab/g2gtools`).
On HPC clusters this usually means loading the `singularity`/`apptainer` module on the compute nodes.

```bash
git clone https://github.com/Itolstoganov/adna_bench.git adna_bench && cd adna_bench
./install.sh            # submodules, adna-eval conda env
conda activate adna-eval
```

## Aligner accuracy benchmark (simx)

1. Adjust datasets, aligner versions and parameters in `simx/config/config.yaml`
   (`simx/config/config_test.yaml` is a small fruitfly config used by `simx/run.py --test`;
   its outputs go to `simx/test/` via the `output_dir` key).

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

`meta_bench/run.sh` runs the benchmark in four steps

1. **databases**: download the [KrakenUniq](https://github.com/fbreitwieser/krakenuniq) and [MALT](https://github.com/husonlab/malt) databases (~750 GB) to `data/` at the repo root. Symlink `data/` to a large disk first if needed.
2. **simulate**: simulate the samples with gargammel using the workflow in `meta_bench/workflow/Simulation.snake`. Composition and parameters are set up in `meta_bench/config/benchmark.yaml`. Reference genomes are in `meta_bench/config/accessions.tsv`. Reads and per-read ground truth will be created in `meta_bench/datasets/sim/{sample}/`.
3. **ameta**: run aMeta ([submodule version](https://github.com/Itolstoganov/aMeta/tree/bowtie2-profiler)) in `meta_bench/runs/sim_dataset/` using configs from `meta_bench/config/ameta.yaml`.
4. **summary**: score the authentication results against the ground truth, and summarise runtime and memory. Results are written to `meta_bench/runs/sim_dataset/results/`. `meta_summary/` holds the ground-truth heatmap and metrics, `meta_summary/bench/` the runtime/memory summary.

```bash
meta_bench/run_full_benchmark.sh               # all four steps
meta_bench/run_full_benchmark.sh 3 4           # only aMeta + summary
meta_bench/run_full_benchmark.sh ameta -- -n   # args after -- go to snakemake
CORES=40 RAM_GB=156 meta_bench/run_full_benchmark.sh
```

**Resources:** The metagenomic benchmark requires ~1.1 TB disk (for krakenuniq and MALT databases), as well as ~160 GB RAM. With 20 cores, MALT takes ~15 min per sample (most of it loading the index); the aligner + [ngsLCA](https://github.com/miwipe/ngsLCA) profilers take a few minutes per sample.

`meta_bench/run.sh --test` runs all four steps on 2 small samples, with databases built from the 6 test genomes (`meta_bench/config/test/`), to check the setup before the full run. It needs ~13 GB disk and ~8 GB RAM.

## Figures 

Output files used in the figures of the [Hybrid RY-mers increase aDNA alignment sensitivity](https://www.biorxiv.org/content/10.64898/2026.11.05.764514v2) paper, produced with default simx and meta_bench configs. Paths are relative to the repository root.

| Figure | Panel | Output file |
|---|---|---|
| 2 | a | `simx/result/summary_plots/outcome.png` |
| 2 | b | `simx/result/summary_plots/outcome_deam.png` |
| 2 | c | `simx/result/summary_plots/accuracy_deam_diff.png` |
| 2 | d | `simx/result/summary_plots/k_runtime_nolegend.png` |
| 2 | e | `simx/result/summary_plots/k_memory.png` |
| 3 | a | `meta_summary/auth_jaccard_threshold.pdf` |
| 3 | b | `meta_summary/bench/bench_profiling_runtime.pdf` |
| 3 | c | `meta_summary/bench/bench_index_runtime.pdf` |
| 3 | d | `meta_summary/bench/bench_peak_memory.pdf` |
| S1 | | `simx/result/summary_plots/combined_k_accuracy.png` |
| S2 | | `simx/result/summary_plots/combined_k_runtime.png` |
| S3 | | `meta_summary/auth_ani_scatter.ancient.pdf` |
| S4 | | `meta_summary/auth_criteria_table.ancient.pdf` |
