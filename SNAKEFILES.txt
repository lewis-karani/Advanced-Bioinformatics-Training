# Define the input fastq files
SAMPLES = ["barcode01","barcode02","barcode03","barcode04","barcode05","barcode06","barcode07","barcode08","barcode09","barcode10","barcode11","barcode12","barcode13","barcode14","barcode15","barcode16","barcode17","barcode18","barcode19","barcode20","barcode21","barcode22","barcode23", "barcode24"] 

rule all:
    input:
        "msa/alignment.fasta",
        "iqtree/alignment.treefile"

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
        trimmomatic SE -threads 23 -phred33 {input} {output} SLIDINGWINDOW:50:10 MINLEN:100
        """

# Step 3: Reference-based Alignment with Minimap2
rule alignment:
    input:
        "data/{sample}_trimmed_filtered.fastq",
        "reference.fasta"
    output:
        "alignment/{sample}_aligned.bam"
    shell:
        """
        minimap2 -ax map-ont {input[1]} {input[0]} | samtools view -bS | samtools sort -o {output}
        samtools index {output}
        """

# Step 4: ivar trimming of primers
rule primertrim:
    input:
        "alignment/{sample}_aligned.bam",
        "ARTIC-V3.bed"
    output:
        "alignment/{sample}_primertrim.bam"
    shell:
        """
        ivar trim -e -i {input[0]} -b {input[1]} | samtools sort -o {output}
        samtools index {output}
        """

# Step 5: Medaka consensus and stitching
rule medaka_consensus:
    input:
        "alignment/{sample}_primertrim.bam"
    output:
        "medaka_output2/{sample}/{sample}.fasta"
    params:
        medaka_dir="medaka_output2/{sample}",
        medaka_hdf="medaka_output2/{sample}/{sample}.hdf"
    shell:
        """
        mkdir -p {params.medaka_dir}
        medaka consensus {input} {params.medaka_hdf} --model r941_min_high_g360 --batch 200 --threads 2
        medaka stitch {params.medaka_hdf} reference.fasta {output}
        """

# Step 6: Update FASTA headers to file names
rule update_fasta_headers:
    input:
        "medaka_output2/{sample}/{sample}.fasta"
    output:
        "medaka_output2/{sample}/{sample}_updated.fasta"
    shell:
        """
        sed -i '1s/.*/>{wildcards.sample}/' {input}
        mv {input} {output}
        """

# Step 7: Concatenate all samples
rule concatenate_consensus:
    input:
        expand("medaka_output2/{sample}/{sample}_updated.fasta" , sample=SAMPLES)
    output:
        "medaka_output2/consensus_concat.fasta"
    shell:
        """
        cat {input} > {output}
        """

# Step 8: Multiple Sequence Alignment with MAFFT
rule multiple_sequence_alignment:
    input:
        "medaka_output2/consensus_concat.fasta"
    output:
        "msa/alignment.fasta"
    shell:
        """
        mafft --auto --reorder {input} > {output}
        """

# Step 9: Build and visualize a phylogenetic tree with IQ-TREE
rule iqtree:
    input:
        "msa/alignment.fasta"
    output:
        "iqtree/alignment.treefile"
    shell:
        """
        mkdir -p iqtree
        iqtree2 -s {input} -m GTR+G -nt AUTO -bb 1000 -pre iqtree/alignment
        """
