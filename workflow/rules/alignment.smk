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

# localrules:
#     result_csv, concat_csvs, concat_genome_csvs, plot_ends, plot_genomes, download_xmapper

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
rule bwa_aln:
    output:
        bam="runs/bwaaln/{genome}-{read_length}/pe.bam",
        sai1=temp("runs/bwaaln/{genome}-{read_length}/1.sai"),
        sai2=temp("runs/bwaaln/{genome}-{read_length}/2.sai")
    input:
        fasta="genomes/{genome}.fa",
        index="genomes/{genome}.fa.bwt",
        r1_fastq="datasets/{genome}/{read_length}/fastp/1.fq.gz",
        r2_fastq="datasets/{genome}/{read_length}/fastp/2.fq.gz"
    params:
        opt="-n 0.01 -o 2 -l 1024"
    threads: 8
    log:
        "runs/bwaaln/{genome}-{read_length}/pe.bam.log"
    shell:
        # "cat {input} > /dev/null; "
        "\n /usr/bin/time -v bwa aln -t {threads} {params.opt} {input.fasta} {input.r1_fastq} > {output.sai1} 2> {log}.tmp"
        "\n /usr/bin/time -v bwa aln -t {threads} {params.opt} {input.fasta} {input.r2_fastq} > {output.sai2} 2> {log}.tmp"
        "\n bwa sampe {input.fasta} {output.sai1} {output.sai2} {input.r1_fastq} {input.r2_fastq}"
        " | grep -v '^@PG'"
        " | samtools view --no-PG -o {output.bam}.tmp.bam -"
        "\n mv -v {output.bam}.tmp.bam {output.bam}"
        "\n mv -v {log}.tmp {log}"


# Map reads with strobealign

rule compile_strobealign:
    output: "bin/strobealign/{commit}"
    threads: 4
    shell:
        """
        wget https://github.com/ksahlin/strobealign/archive/{wildcards.commit}.zip
        unzip -d strobealign-{wildcards.commit} {wildcards.commit}.zip
        rm {wildcards.commit}.zip
        cd strobealign-{wildcards.commit}/*
        cmake -B build -DCMAKE_C_FLAGS="-march=native" -DCMAKE_CXX_FLAGS="-march=native"
        make -j {threads} -C build strobealign
        mv build/strobealign ../../{output}
        cd ../..
        rm -rf strobealign-{wildcards.commit}
        """

rule strobealign_paired_end:
    output:
        bam="runs/strobealign-{program}/{genome}-{read_length}/pe.bam"
    input:
        binary=lambda wildcards: VERSIONS[wildcards.program]["binary"],
        fasta="genomes/{genome}.fa",
        r1_fastq="datasets/{genome}/{read_length}/fastp/1.fq.gz",
        r2_fastq="datasets/{genome}/{read_length}/fastp/2.fq.gz",
    params:
        extra_args=lambda wildcards: VERSIONS[wildcards.program]["arguments"]
    threads: 8
    log:
        "runs/strobealign-{program}/{genome}-{read_length}/pe.bam.log"
    shell:
        "cat {input} > /dev/null; "
        "/usr/bin/time -v {input.binary}{params.extra_args} -v -t {threads} {input.fasta} {input.r1_fastq} {input.r2_fastq} 2> {log}.tmp"
        " | grep -v '^@PG'"
        " | samtools view --no-PG -o {output.bam}.tmp.bam -"
        "\n mv -v {output.bam}.tmp.bam {output.bam}"
        "\n mv -v {log}.tmp {log}"


# rule strobealign_single_end:
#     output:
#         bam=temp("runs/strobealign-{program}/{dataset}-{genome}-{read_length}/se.bam")
#     input:
#         binary=lambda wildcards: VERSIONS[wildcards.program]["binary"],
#         fasta="genomes/{genome}.fa",
#         r1_fastq="datasets/{dataset}/{genome}-{read_length}/1.fastq.gz",
#     params:
#         extra_args=lambda wildcards: VERSIONS[wildcards.program]["arguments"]
#     threads: 8
#     log:
#         "runs/strobealign-{program}/{dataset}-{genome}-{read_length}/se.bam.log"
#     shell:
#         "cat {input} > /dev/null; "
#         "/usr/bin/time -v {input.binary}{params.extra_args} -v -t {threads} {input.fasta} {input.r1_fastq} 2> {log}.tmp"
#         " | grep -v '^@PG'"
#         " | samtools view --no-PG -o {output.bam}.tmp.bam -"
#         "\n mv -v {output.bam}.tmp.bam {output.bam}"
#         "\n mv -v {log}.tmp {log}"


rule resources:
    output:
        csv="csv/{prog}/{genome}-{read_length}/{ends}.{bampaf}.resources.csv"
    input:
        log="runs/{prog}/{genome}-{read_length}/{ends}.{bampaf}.log",
    params:
        typ=lambda wildcards: "map" if wildcards.bampaf == "paf" else "align"
    shell:
        """
        echo -n {wildcards.prog},{params.typ},{wildcards.genome},{wildcards.ends},{wildcards.read_length}, > {output.csv}.tmp
        user_time=$(awk '/User time/ {{print $NF}}' {input.log})
        memory=$(awk '/Maximum resident/ {{print $NF/1048576}}' {input.log})
        indexing_time=$(awk '/Total time indexing:/ {{print $4}}' {input.log})
        mapping_time=$(awk '/Total time mapping:/ {{print $4}}' {input.log})
        if test -z "${{mapping_time}}"; then mapping_time=$(sed -n -e 's|.* Real time: \\([1-9][0-9.]*\\) sec.*|\\1|p' {input.log}); fi
        if test -z "${{mapping_time}}"; then mapping_time=$(awk '/^Processing query/ && $3 >= 10000 && t == 0 {{t=substr($5, 1, length($5)-1);}}; /^Done in/ {{total=substr($3, 1, length($3)-2)}}; END{{print total-t}}' {input.log}); fi
        echo ,${{mapping_time}},${{memory}} >> {output}.tmp
        mv {output.csv}.tmp {output.csv}
        """
        

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
