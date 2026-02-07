# Alignment rules are largely taken from https://github.com/NBISweden/strobealign-evaluation

# MODES = config.get("modes", ("align", "map"))
# ENDS = config.get("ends", ("se", "pe"))
# DATASETS = config["datasets"]
# GENOMES = config["genomes"]
# READ_LENGTHS = config["read-lengths"]

# VERSIONS = {}
# for version in config.get("versions", []):
#     key = version["commit"]
#     if arg := version.get("arguments", ""):
#         key += "-" + arg.replace(" ", "").replace("-", "").replace("=", "")
#     VERSIONS[key] = {
#         "commit": version["commit"],
#         "arguments": f" {arg}" if arg else "",
#         "name": version["name"],
#         "color": version["color"],
#         "binary": f"bin/strobealign/{version['commit']}",
#     }

# PROGRAMS = expand("strobealign-{key}", key=VERSIONS.keys())
# if config["programs"]:
#     PROGRAMS.extend(config["programs"])
# PLOTS = ("accuracy", "aligned", "time", "memory")
# For reads shorter than this, the evaluation includes paired-end data
LONG_READ_LENGTH = 1000

import itertools

localrules:
    safari_install

# rule final:
#     input:
#         expand("plots/ends-{ends}-{plot}.pdf", ends=ENDS, plot=PLOTS),
#         expand("plots/genome-{genome}-{plot}.pdf", genome=GENOMES, plot=PLOTS),
#         "result.csv"


# Map reads with bwa aln

rule bwa_index:
    output:
        "{genome}.fa.amb",
        "{genome}.fa.ann",
        "{genome}.fa.bwt",
        "{genome}.fa.pac",
        "{genome}.fa.sa"
    input: "{genome}.fa"
    shell:
        "bwa index {input}"

# Fix runtime measurements
rule bwa_aln_paired_end:
    output:
        bam="runs/bwaaln/{dataset_id}/pe.bam",
        sai1=temp("runs/bwaaln/{dataset_id}/1.sai"),
        sai2=temp("runs/bwaaln/{dataset_id}/2.sai")
    input:
        fasta="datasets/{dataset_id}/ref.fa",
        index="datasets/{dataset_id}/ref.fa.bwt",
        r1_fastq="datasets/{dataset_id}/fastp/1.fq.gz",
        r2_fastq="datasets/{dataset_id}/fastp/2.fq.gz"
    params:
        opt="-n 0.01 -o 2 -l 1024"
    threads: 8
    log:
        "runs/bwaaln/{dataset_id}/pe.bam.log"
    shell:
        # "cat {input} > /dev/null; "
        "\n /usr/bin/time -v bwa aln -t {threads} {params.opt} {input.fasta} {input.r1_fastq} > {output.sai1} 2> {log}.tmp"
        "\n /usr/bin/time -v bwa aln -t {threads} {params.opt} {input.fasta} {input.r2_fastq} > {output.sai2} 2> {log}.tmp"
        "\n bwa sampe {input.fasta} {output.sai1} {output.sai2} {input.r1_fastq} {input.r2_fastq}"
        " | grep -v '^@PG'"
        " | samtools view --no-PG -o {output.bam}.tmp.bam -"
        "\n mv -v {output.bam}.tmp.bam {output.bam}"
        "\n mv -v {log}.tmp {log}"


rule bwa_aln_single_end:
    output:
        bam="runs/bwaaln/{dataset_id}/se.bam",
        sai1=temp("runs/bwaaln/{dataset_id}/se.sai")
        # sai2=temp("runs/bwaaln/{dataset_id}/2.sai")
    input:
        fasta="datasets/{dataset_id}/ref.fa",
        index="datasets/{dataset_id}/ref.fa.bwt",
        r1_fastq="datasets/{dataset_id}/fastp/1.fq.gz"
    params:
        opt="-n 0.01 -o 2 -l 1024"
    threads: 8
    log:
        "runs/bwaaln/{dataset_id}/se.bam.log"
    shell:
        # "cat {input} > /dev/null; "
        "\n /usr/bin/time -v bwa aln -t {threads} {params.opt} {input.fasta} {input.r1_fastq} > {output.sai1} 2> {log}.tmp"
        "\n bwa samse {input.fasta} {output.sai1} {input.r1_fastq} "
        " | grep -v '^@PG'"
        " | samtools view --no-PG -o {output.bam}.tmp.bam -"
        "\n mv -v {output.bam}.tmp.bam {output.bam}"
        "\n mv -v {log}.tmp {log}"


rule bwa_mem_paired_end:
    output:
        bam="runs/bwamem/{dataset_id}/pe.bam"
    input:
        fasta="datasets/{dataset_id}/ref.fa",
        index="datasets/{dataset_id}/ref.fa.bwt",
        r1_fastq="datasets/{dataset_id}/fastp/1.fq.gz",
        r2_fastq="datasets/{dataset_id}/fastp/2.fq.gz",
    threads: 8
    log:
        "runs/bwamem/{dataset_id}/pe.bam.log"
    shell:
        "cat {input} > /dev/null; "
        "/usr/bin/time -v bwa mem -t {threads} {input.fasta} {input.r1_fastq} {input.r2_fastq} 2> {log}.tmp"
        " | grep -v '^@PG'"
        " | samtools view --no-PG -o {output.bam}.tmp.bam -"
        "\n mv -v {output.bam}.tmp.bam {output.bam}"
        "\n mv -v {log}.tmp {log}"


rule bwa_mem_single_end:
    output:
        bam="runs/bwamem/{dataset_id}/se.bam"
    input:
        fasta="datasets/{dataset_id}/ref.fa",
        index="datasets/{dataset_id}/ref.fa.bwt",
        r1_fastq="datasets/{dataset_id}/fastp/1.fq.gz"
    threads: 8
    log:
        "runs/bwamem/{dataset_id}/se.bam.log"
    shell:
        "cat {input} > /dev/null; "
        "/usr/bin/time -v bwa mem -t {threads} {input.fasta} {input.r1_fastq} 2> {log}.tmp"
        " | grep -v '^@PG'"
        " | samtools view --no-PG -o {output.bam}.tmp.bam -"
        "\n mv -v {output.bam}.tmp.bam {output.bam}"
        "\n mv -v {log}.tmp {log}"


# Map reads with strobealign

# rule compile_strobealign:
#     output: "bin/strobealign/{commit}"
#     threads: 4
#     shell:
#         """
#         wget https://github.com/ksahlin/strobealign/archive/{wildcards.commit}.zip
#         unzip -d strobealign-{wildcards.commit} {wildcards.commit}.zip
#         rm {wildcards.commit}.zip
#         cd strobealign-{wildcards.commit}/*
#         cmake -B build -DCMAKE_C_FLAGS="-march=native" -DCMAKE_CXX_FLAGS="-march=native"
#         make -j {threads} -C build strobealign
#         mv build/strobealign ../../{output}
#         cd ../..
#         rm -rf strobealign-{wildcards.commit}
#         """

# todo release/debug option?
rule compile_strobealign:
    output: "bin/strobealign/{commit}"
    threads: 4
    shell:
        """
        wget https://github.com/ksahlin/strobealign/archive/{wildcards.commit}.zip
        unzip -d strobealign-{wildcards.commit} {wildcards.commit}.zip
        rm {wildcards.commit}.zip
        cd strobealign-{wildcards.commit}/*
        cargo build -j 4
        mv target/debug/strobealign ../../{output}
        cd ../..
        rm -rf strobealign-{wildcards.commit}
        """

rule strobealign_paired_end:
    output:
        bam="runs/strobealign-{program}/{dataset_id}/pe.bam"
    input:
        binary=lambda wildcards: VERSIONS[wildcards.program]["binary"],
        fasta="datasets/{dataset_id}/ref.fa",
        r1_fastq="datasets/{dataset_id}/fastp/1.fq.gz",
        r2_fastq="datasets/{dataset_id}/fastp/2.fq.gz",
    params:
        extra_args=lambda wildcards: VERSIONS[wildcards.program]["arguments"]
    threads: 8
    log:
        "runs/strobealign-{program}/{dataset_id}/pe.bam.log"
    shell:
        "cat {input} > /dev/null; "
        "/usr/bin/time -v {input.binary}{params.extra_args} -v -t {threads} {input.fasta} {input.r1_fastq} {input.r2_fastq} 2> {log}.tmp"
        " | grep -v '^@PG'"
        " | samtools view --no-PG -o {output.bam}.tmp.bam -"
        "\n mv -v {output.bam}.tmp.bam {output.bam}"
        "\n mv -v {log}.tmp {log}"


rule strobealign_single_end:
    output:
        bam="runs/strobealign-{program}/{dataset_id}/se.bam"
    input:
        binary=lambda wildcards: VERSIONS[wildcards.program]["binary"],
        fasta="datasets/{dataset_id}/ref.fa",
        r1_fastq="datasets/{dataset_id}/fastp/1.fq.gz"
    params:
        extra_args=lambda wildcards: VERSIONS[wildcards.program]["arguments"]
    threads: 8
    log:
        "runs/strobealign-{program}/{dataset_id}/se.bam.log"
    shell:
        "cat {input} > /dev/null; "
        "/usr/bin/time -v {input.binary}{params.extra_args} -v -t {threads} {input.fasta} {input.r1_fastq} 2> {log}.tmp"
        " | grep -v '^@PG'"
        " | samtools view --no-PG -o {output.bam}.tmp.bam -"
        "\n mv -v {output.bam}.tmp.bam {output.bam}"
        "\n mv -v {log}.tmp {log}"


rule safari_install:
    output:
        bin="bin/safari/vg"
    shell:
        """
        wget https://github.com/grenaud/SAFARI/releases/download/v1.0.1/vg
        mv vg {output.bin}
        chmod 755 {output.bin}
        """


rule safari_index_linear:
    #todo are all these really necessary
    #todo add runtime/mem?
    input:
        safari="bin/safari/vg",
        fasta="datasets/{dataset_id}/ref.fa"
    output:
        vg="datasets/{dataset_id}/safari/ref.vg",
        og="datasets/{dataset_id}/safari/ref.og",
        gfa="datasets/{dataset_id}/safari/ref.gfa",
        dist="datasets/{dataset_id}/safari/ref.dist",
        xg="datasets/{dataset_id}/safari/ref.xg",
        gbz="datasets/{dataset_id}/safari/ref.gbz",
        gg="datasets/{dataset_id}/safari/ref.gg",
        gbwt="datasets/{dataset_id}/safari/ref.gbwt",
        min="datasets/{dataset_id}/safari/ref.min",
        ry="datasets/{dataset_id}/safari/ref.ry"
    threads:
        16
    shell:
        """
        {input.safari} construct -p -r {input.fasta} -t {threads} > {output.vg}
        {input.safari} convert -t {threads} -o {output.vg} > {output.og}
        {input.safari} view {output.og} > {output.gfa}
        {input.safari} index -t {threads} -j {output.dist} {output.og}
        {input.safari} index -t {threads} -x {output.xg} {output.og}
        {input.safari} gbwt --num-threads {threads} -g {output.gbz} --gbz-format -G {output.gfa}
        {input.safari} gbwt --num-threads {threads} -o {output.gbwt} -g {output.gg} -Z {output.gbz}
        {input.safari} minimizer -t {threads} -o {output.min} -g {output.gbwt} -d {output.dist} {output.og}
        {input.safari} rymer -t {threads} -o {output.ry} -g {output.gbwt} -d {output.dist} {output.og}
        """

rule safari_single_end:
    # todo use more suitable damage profiles
    input:
        fasta="datasets/{dataset_id}/ref.fa",
        deam_3p="safari_data/jeong_2020_cell_SHG003_MT_3p.prof",
        deam_5p="safari_data/jeong_2020_cell_SHG003_MT_5p.prof",
        r1_fastq="datasets/{dataset_id}/fastp/1.fq.gz",
        min="datasets/{dataset_id}/safari/ref.min",
        dist="datasets/{dataset_id}/safari/ref.dist",
        ry="datasets/{dataset_id}/safari/ref.ry",
        gbz="datasets/{dataset_id}/safari/ref.gbz",
        safari="bin/safari/vg"
    output:
        bam="runs/safari/{dataset_id}/se.bam"
    log:
        "runs/safari/{dataset_id}/se.bam.log"
    params:
        j=0.5
    threads:
        8
    shell:
        "cat {input} > /dev/null; "
        "/usr/bin/time -v {input.safari} safari -t {threads} -f {input.r1_fastq} -j {params.j}"
        " --deam-3p {input.deam_3p} --deam-5p {input.deam_5p}" 
        " -m {input.min} -d {input.dist} -q {input.ry} -Z {input.gbz} -o bam 2> {log}.tmp"
        # " | grep -v '^@PG'"
        # " | samtools view -o {output.bam}.tmp.nomd.bam "
        " | samtools calmd -b - {input.fasta} --no-PG > {output.bam}.tmp.bam"
        "\n mv -v {output.bam}.tmp.bam {output.bam}"
        "\n mv -v {log}.tmp {log}"


rule safari_paired_end:
    # todo use more suitable damage profiles
    input:
        fasta="datasets/{dataset_id}/ref.fa",
        deam_3p="safari_data/jeong_2020_cell_SHG003_MT_3p.prof",
        deam_5p="safari_data/jeong_2020_cell_SHG003_MT_5p.prof",
        r1_fastq="datasets/{dataset_id}/fastp/1.fq.gz",
        r2_fastq="datasets/{dataset_id}/fastp/2.fq.gz",
        min="datasets/{dataset_id}/safari/ref.min",
        dist="datasets/{dataset_id}/safari/ref.dist",
        ry="datasets/{dataset_id}/safari/ref.ry",
        gbz="datasets/{dataset_id}/safari/ref.gbz",
        safari="bin/safari/vg"
    output:
        bam="runs/safari/{dataset_id}/pe.bam"
    log:
        "runs/safari/{dataset_id}/pe.bam.log"
    params:
        j=0.5
    threads:
        8
    shell:
        "cat {input} > /dev/null; "
        "/usr/bin/time -v {input.safari} safari -t {threads} -f {input.r1_fastq} -f {input.r2_fastq} -j {params.j}"
        " --deam-3p {input.deam_3p} --deam-5p {input.deam_5p}" 
        " -m {input.min} -d {input.dist} -q {input.ry} -Z {input.gbz} -o bam 2> {log}.tmp"
        # " | grep -v '^@PG'"
        # " | samtools view -o {output.bam}.tmp.nomd.bam "
        " | samtools calmd -b - {input.fasta} --no-PG > {output.bam}.tmp.bam"
        "\n mv -v {output.bam}.tmp.bam {output.bam}"
        "\n mv -v {log}.tmp {log}"

# def csvs(wildcards):
#     files = []
#     if "align" in MODES:
#         for prog, dataset, read_length in itertools.product(PROGRAMS, DATASETS, READ_LENGTHS):
#             if read_length >= LONG_READ_LENGTH and wildcards.ends == "pe":
#                 continue
#             if prog == "xmapper" and read_length >= LONG_READ_LENGTH:
#                 continue
#             if prog == "xmapper" and wildcards.genome == "rye":
#                 # X-Mapper crashes on rye
#                 continue
#             files.append(f"csv/{prog}/{dataset}-{wildcards.genome}-{read_length}/{wildcards.ends}.bam.csv")
#     if "map" in MODES:
#         if wildcards.ends == "se":
#                read_lengths = READ_LENGTHS
#         else:
#             read_lengths = [rl for rl in READ_LENGTHS if rl < LONG_READ_LENGTH]
#         files += expand(
#             f"csv/{{prog}}/{{dataset}}-{wildcards.genome}-{{read_length}}/{wildcards.ends}.paf.csv",
#             prog=(prog for prog in PROGRAMS if prog not in {"bwamem", "xmapper"}), dataset=DATASETS, read_length=read_lengths
#         )
#     return files


# rule concat_csvs:
#     output: "csv/{genome,([^/]*)}-{ends,(se|pe)}.csv"
#     input: csvs
#     shell:
#         "( echo tool,mode,dataset,genome,ends,read_length,read_count,aligned,overaligned,accuracy,jaccuracy,saccuracy,time,memory ; cat {input} ) > {output}"


# rule concat_genome_csvs:
#     output: "result.csv"
#     input: expand("csv/{genome}-{ends}.csv", ends=ENDS, genome=GENOMES)
#     shell:
#         "( head -n 1 {input[0]} ; for f in {input}; do sed 1d $f; done ) > {output}"

# rule plot_ends:
#     output: expand("plots/ends-{end}-{plot}.pdf", end=ENDS, plot=PLOTS)
#     input:
#         csv="result.csv",
#         config="config.yaml",
#     shell:
#         "python plots.py -c {input.config} {input.csv} plots/"

# rule plot_genomes:
#     output: expand("plots/genome-{genome}-{plot}.pdf", genome=GENOMES, plot=PLOTS)
#     input:
#         csv="result.csv",
#         config="config.yaml",
#     shell:
#         "python plots.py --genome -c {input.config} {input.csv} plots/"
