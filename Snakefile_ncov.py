# Snakefile for SARS-CoV-2 nanopore reads analysis to generate a phylogenetic tree

# Define the input fastq files
samples = ["*.fastq.gz"]
fastq_files = expand("data/barcode*/{sample}", sample=samples)

rule all:
    input:
        "results/phylogenetic_tree.nwk"

# Step 1: Adapter Trimming with Porechop
rule porechop:
    input:
        "data/{sample}.fastq"
    output:
        "data/{sample}_trimmed.fastq"
    shell:
        """
        porechop -i {input} -o {output}
        """

# Step 2: Quality Trimming with Trimmomatic
rule trimmomatic:
    input:
        "data/{sample}_trimmed.fastq"
    output:
        "data/{sample}_trimmed_filtered.fastq"
  shell:
        """
        trimmomatic SE {input} {output} SLIDINGWINDOW:50:10 MINLEN:100
        """

# Step 3: Reference-based Alignment with Minimap2
rule alignment:
    input:
        "data/{sample}_trimmed_filtered.fastq",
        "reference/reference.fasta"
    output:
        "alignment/{sample}_aligned.bam"
   shell:
        """
        minimap2 -ax map-ont {input[1]} {input[0]} | samtools view -bS | samtools sort -o {output}
        samtools index {output}
        """

# Step 4: Variant Calling and Consensus Generation with iVar
rule ivar:
    input:
        bam="alignment/{sample}_aligned.bam",
        bai="alignment/{sample}_aligned.bam.bai",
        ref="reference/reference.fasta"
    output:
        "variants/{sample}.vcf",
        "consensus/{sample}_consensus.fa"
    shell:
        """
        samtools mpileup -A -d 0 -Q 0 -u -v -f {input.ref} {input.bam} | bcftools call -cv -Oz -o {output[0]}
        bcftools index {output[0]}
        ivar consensus -p {output[1]} -q 20 -t 0.5 -m 10 -n N -i {input.ref} -b {input.bam} -v {output[0]}
        """

# Step 5: Polishing with Medaka
rule medaka:
    input:
        bam="alignment/{sample}_aligned.bam",
        ref="reference/sars_cov_2_reference.fasta",
        consensus="consensus/{sample}_consensus.fa"
    output:
        "polished/{sample}_polished.fasta"
    shell:
        """
        medaka_consensus -i {input.bam} -d {input.consensus} -o polished/{wildcards.sample}
        cp polished/{wildcards.sample}/consensus.fasta {output}
        """

# Step 6: Multiple Sequence Alignment with MAFFT
rule multiple_sequence_alignment:
    input:
        expand("polished/{sample}_polished.fasta", sample=samples)
    output:
        "msa/alignment.fasta"
   shell:
        """
        mafft --auto --reorder {input} > {output}
        """

# Step 7: Phylogenetic Tree Construction with IQ-TREE
rule phylogenetic_tree:
    input:
        "msa/alignment.fasta"
    output:
        "results/phylogenetic_tree.nwk"
   shell:
        """
        iqtree -s {input} -nt AUTO -m GTR+G -bb 1000 -alrt 1000 -pre results/phylogenetic_tree
        cp results/phylogenetic_tree.treefile {output}
        """
