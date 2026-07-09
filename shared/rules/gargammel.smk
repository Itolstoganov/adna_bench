# Shared gargammel simulation rules, reused by sim_bench (workflow/) and
# meta_bench (meta_bench/workflow/). Both run from the repo root.
#
# Contract: before `include:`-ing this file, the including workflow must define a
# dict `GARGAMMEL_REGISTRY` mapping an output prefix ("simdir") to a params dict:
#
#   GARGAMMEL_REGISTRY["datasets/foo"] = {
#       "endo": ["genomes/foo.fa"],        # list of endogenous FASTAs
#       "cont": ["genomes/bar.fa"],        # list of contaminant FASTAs
#       "bact_dir": "gargammel_data/k14/fasta",   # sim_bench: symlink dir as bact/
#       # --- OR, for explicit composition (meta_bench): ---
#       "bact_list": "…/list", "bact_genome_dir": "meta_bench/genomes",
#       "num": 500000, "comp": "0.2,0.0,0.8",
#       "gargammel_args": "-l 40 -damage 0.03,0.4,0.01,0.3",
#       "extra_inputs": [ … ],              # extra files the sim rule depends on
#   }
#
# (Ground-truth parsing is done by each workflow's own thin rule, not here.)
#
# gargammel always writes both simulated_s1 and simulated_s2 (single-end callers
# simply consume s1). This file provides ONLY the shared `gargammel_simulate`
# rule (the heavy, otherwise-duplicated compute). Ground-truth parsing differs
# per workflow (sim_bench parses fastp-trimmed reads; meta_bench parses the raw
# per-pass reads with taxon labels), so each workflow keeps a thin parse rule
# that calls the shared shared/scripts/parse_gargammel.py.
#   Output per simdir: {simdir}/simulated_s{1,2}.fq.gz

import re

SHARED_DIR = "shared"
GARGAMMEL_PL = config.get("gargammel_pl", "gargammel/gargammel.pl")

# Constrain the generic {simdir} wildcard to registered prefixes so this rule
# never competes ambiguously with other rules.
if GARGAMMEL_REGISTRY:
    wildcard_constraints:
        simdir="(" + "|".join(re.escape(k) for k in GARGAMMEL_REGISTRY) + ")",


def _greg(wildcards):
    return GARGAMMEL_REGISTRY[wildcards.simdir]


def gargammel_sim_inputs(wildcards):
    r = _greg(wildcards)
    inputs = [GARGAMMEL_PL]
    inputs += list(r.get("endo", []))
    inputs += list(r.get("cont", []))
    if r.get("bact_list"):
        inputs.append(r["bact_list"])
    inputs += list(r.get("extra_inputs", []))
    return inputs


rule gargammel_simulate:
    input:
        gargammel_sim_inputs,
    output:
        s1="{simdir}/simulated_s1.fq.gz",
        s2="{simdir}/simulated_s2.fq.gz",
    params:
        run=f"{SHARED_DIR}/scripts/run_gargammel.py",
        gargammel=GARGAMMEL_PL,
        num=lambda wc: _greg(wc)["num"],
        comp=lambda wc: _greg(wc)["comp"],
        gargammel_args=lambda wc: _greg(wc).get("gargammel_args", ""),
        outprefix=lambda wc: f"{wc.simdir}/simulated",
        endo=lambda wc: " ".join(f"--endo {p}" for p in _greg(wc).get("endo", [])),
        cont=lambda wc: " ".join(f"--cont {p}" for p in _greg(wc).get("cont", [])),
        bact=lambda wc: (
            f"--bact-dir {_greg(wc)['bact_dir']}"
            if _greg(wc).get("bact_dir")
            else f"--bact-list {_greg(wc)['bact_list']} "
                 f"--bact-genome-dir {_greg(wc)['bact_genome_dir']}"
        ),
    shell:
        r"""
        python {params.run} \
            --workdir {wildcards.simdir} \
            --out {params.outprefix} \
            --gargammel {params.gargammel} \
            --num {params.num} --comp {params.comp} \
            {params.endo} {params.cont} {params.bact} \
            --gargammel-args {params.gargammel_args:q}
        """
