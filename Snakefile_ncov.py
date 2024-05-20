# Snakefile

import os

RAW_READS_DIR = "raw_reads"
ADAPTER_REMOVED_DIR = "adapter_removed"
TRIMMED_READS_DIR = "trimmed_reads"
READMAPPING_DIR = "readmapping"
REFERENCE = "../Reference/sarscov2-Wu1.fasta"
BED_FILE = "../artic_v3/ARTIC-V3.bed"
MEDAKA_OUTPUT_DIR = "./medaka_output2"
THREADS = 23

rule all:
    input:
        expand(READMAPPING_DIR + "/{sample}_primertrim_sorted.bam", sample=[f"sample{i:02d}" for i in range(1, 25)]),
        expand(MEDAKA_OUTPUT_DIR + "/sample{sample:02d}/sample{sample:02d}.fasta", sample=range(1, 25))

rule porechop:
    input:
        RAW_READS_DIR + "/{sample}.fastq.gz"
    output:
        ADAPTER_REMOVED_DIR + "/{sample}_adapter_rm.fastq.gz"
    shell:
        "porechop -i {input} -o {output}"

rule trimmomatic:
    input:
        ADAPTER_REMOVED_DIR + "/{sample}_adapter_rm.fastq.gz"
    output:
        TRIMMED_READS_DIR + "/{sample}_trimmed.fastq.gz"
    shell:
        "trimmomatic SE -threads {THREADS} -phred33 {input} {output} SLIDINGWINDOW:50:10 MINLEN:100"

rule unpigz:
    input:
        TRIMMED_READS_DIR + "/{sample}_trimmed.fastq.gz"
    output:
        TRIMMED_READS_DIR + "/{sample}_trimmed.fastq"
    shell:
        "unpigz {input}"

rule minimap2:
    input:
        TRIMMED_READS_DIR + "/{sample}_trimmed.fastq"
    output:
        READMAPPING_DIR + "/{sample}_aln.bam"
    shell:
        "minimap2 -ax map-ont {REFERENCE} {input} > {output}"

rule samtools_sort:
    input:
        READMAPPING_DIR + "/{sample}_aln.bam"
    output:
        READMAPPING_DIR + "/{sample}_sorted.bam"
    shell:
        "samtools sort {input} > {output}"

rule ivar_trim:
    input:
        READMAPPING_DIR + "/{sample}_sorted.bam"
    output:
        READMAPPING_DIR + "/{sample}_primertrim.bam"
    shell:
        "ivar trim -e -i {input} -b {BED_FILE} -p {output.replace('.bam', '')}"

rule samtools_sort_trimmed:
    input:
        READMAPPING_DIR + "/{sample}_primertrim.bam"
    output:
        READMAPPING_DIR + "/{sample}_primertrim_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"

rule medaka_consensus:
    input:
        READMAPPING_DIR + "/{sample}_primertrim_sorted.bam"
    output:
        MEDAKA_OUTPUT_DIR + "/{sample}/{sample}.fasta"
    params:
        hdf="{MEDAKA_OUTPUT_DIR}/{sample}/{sample}.hdf",
        output_dir=lambda wildcards: f"{MEDAKA_OUTPUT_DIR}/{wildcards.sample}"
    run:
        shell("mkdir -p {params.output_dir}")
        shell("medaka consensus {input} {params.hdf} --model r941_min_high_g360 --batch 200 --threads 2")
        shell("medaka stitch {params.hdf} {REFERENCE} {output}")
