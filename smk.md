# Snakefile for SARS-CoV-2 nanopore reads analysis to generate a phylogenetic tree

# Define the conda environment for each rule
envs = {
    "basecalling": "envs/basecalling.yml",
    "quality_control": "envs/quality_control.yml",
    "assembly": "envs/assembly.yml",
    "alignment": "envs/alignment.yml",
    "phylogeny": "envs/phylogeny.yml"
}

# Define the input fastq files
samples = ["sample1", "sample2", "sample3"]
fastq_files = expand("data/{sample}.fastq", sample=samples)

rule all:
    input:
        "results/phylogenetic_tree.nwk"

# Step 1: Basecalling
rule basecalling:
    input:
        "data/raw_reads/{sample}.fast5"
    output:
        "data/{sample}.fastq"
    conda:
        envs["basecalling"]
    shell:
        """
        guppy_basecaller -i {input} -s {output} --config dna_r9.4.1_450bps_fast.cfg
        """

# Step 2: Quality Control
rule quality_control:
    input:
        "data/{sample}.fastq"
    output:
        "data/{sample}_filtered.fastq"
    conda:
        envs["quality_control"]
    shell:
        """
        nanofilt -q 7 --length 400 < {input} > {output}
        """

# Step 3: Genome Assembly
rule genome_assembly:
    input:
        "data/{sample}_filtered.fastq"
    output:
        "assembly/{sample}_assembly.fasta"
    conda:
        envs["assembly"]
    shell:
        """
        flye --nano-raw {input} --out-dir assembly/{wildcards.sample} --genome-size 30k
        """

# Step 4: Reference-based Alignment
rule alignment:
    input:
        "assembly/{sample}_assembly.fasta",
        "reference/sars_cov_2_reference.fasta"
    output:
        "alignment/{sample}_aligned.sam"
    conda:
        envs["alignment"]
    shell:
        """
        minimap2 -a {input[1]} {input[0]} > {output}
        """

# Step 5: Consensus Generation
rule consensus:
    input:
        "alignment/{sample}_aligned.sam"
    output:
        "consensus/{sample}_consensus.fasta"
    conda:
        envs["alignment"]
    shell:
        """
        samtools view -Sb {input} | samtools sort -o {wildcards.sample}.sorted.bam
        samtools mpileup -uf reference/sars_cov_2_reference.fasta {wildcards.sample}.sorted.bam | bcftools call -mv -Oz -o {wildcards.sample}.vcf.gz
        bcftools index {wildcards.sample}.vcf.gz
        bcftools consensus -f reference/sars_cov_2_reference.fasta {wildcards.sample}.vcf.gz > {output}
        """

# Step 6: Multiple Sequence Alignment
rule multiple_sequence_alignment:
    input:
        expand("consensus/{sample}_consensus.fasta", sample=samples)
    output:
        "msa/alignment.fasta"
    conda:
        envs["alignment"]
    shell:
        """
        mafft --auto --reorder {input} > {output}
        """

# Step 7: Phylogenetic Tree Construction
rule phylogenetic_tree:
    input:
        "msa/alignment.fasta"
    output:
        "results/phylogenetic_tree.nwk"
    conda:
        envs["phylogeny"]
    shell:
        """
        fasttree -gtr -nt {input} > {output}
        """
