# Reference-genome acquisition and per-sample composition for the metagenomic
# simulation. Included by meta_bench/workflow/Simulation.snake (run from repo root).

MB_SCRIPTS = "meta_bench/workflow/scripts"

PASS_TABLE = {"ancient": config["ancient_table"], "modern": config["modern_table"]}


rule download_genome:
    output:
        fa=GENOMES + "/{species}.fa",
    params:
        accession=lambda wc: SPECIES[wc.species]["accession"],
        script=f"{MB_SCRIPTS}/fetch_genome.py",
    wildcard_constraints:
        species="(" + "|".join(re.escape(s) for s in SPECIES) + ")",
    shell:
        r"""
        python {params.script} --species {wildcards.species} \
            --accession '{params.accession}' --out {output.fa}
        """


rule build_contig2taxon:
    input:
        genomes=expand(GENOMES + "/{species}.fa", species=MICROBIAL_SPECIES),
    output:
        tsv=GENOMES + "/contig2taxon.tsv",
    params:
        script=f"{MB_SCRIPTS}/build_contig2taxon.py",
        accessions=config["accessions"],
    shell:
        r"""
        python {params.script} --genomes {input.genomes} \
            --accessions {params.accessions} --out {output.tsv}
        """


rule build_list:
    output:
        lst=RESULTS + "/{sample}/{phase}/bact.list",
    params:
        script=f"{MB_SCRIPTS}/build_lists.py",
        table=lambda wc: PASS_TABLE[wc.phase],
        sample_num=lambda wc: int(wc.sample.replace("sample", "")),
    wildcard_constraints:
        phase="(ancient|modern)",
        sample=r"sample\d+",
    shell:
        r"""
        python {params.script} --table {params.table} \
            --sample {params.sample_num} --out {output.lst}
        """
