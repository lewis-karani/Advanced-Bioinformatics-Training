cd raw_reads/

for file in *.fastq.gz; do
porechop -i "$file" -o "${file//.fastq.gz/_adapter_rm.fastq.gz}"
done
mv *_adapter_rm.fastq.gz ../adapter_removed/
cd ../adapter_removed
for file in *.fastq.gz; do
trimmomatic SE -threads 23 -phred33 $file ${file//_adapter_rm/_trimmed} SLIDINGWINDOW:50:10 MINLEN:100
done
mv *_trimmed.fastq.gz ../trimmed_reads/
cd ../trimmed_reads
unpigz *.fastq.gz
for file in *trimmed.fastq; do
minimap2 -ax map-ont ../Reference/sarscov2-Wu1.fasta $file > "${file//.fastq/_aln.bam}"
done
mv *aln.bam ../readmapping
pigz *.fastq
cd ../readmapping
for file in *aln.bam; do
samtools sort $file > "${file//.bam/_sorted.bam}"
done
for file in *sorted.bam; do
ivar trim -e -i $file -b ../artic_v3/ARTIC-V3.bed -p "${file//sorted.bam/primertrim.bam}" 
done
for file in *primertrim.bam; do
samtools sort $file -o "${file//primertrim/primertrim_sorted}"
done
for file in readmapping/*primertrim_sorted.bam; do
    if [ -f "$file" ]; then
        filename_without_extension="$(basename -- "$file")"
        filename_without_extension="${filename_without_extension%.*}"
        mkdir -p "./medaka_output2/$filename_without_extension"
        medaka consensus $file medaka_output2/$filename_without_extension/$filename_without_extension.hdf --model r941_min_high_g360 --batch 200 --threads 2
        medaka stitch medaka_output2/$filename_without_extension/$filename_without_extension.hdf sarscov2-Wu1.fasta medaka_output2/$filename_without_extension/$filename_without_extension.fasta
    fi
done