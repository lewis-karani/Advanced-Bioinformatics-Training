for i in {01..24}; do cd barcode$i/ && unpigz * && cat *.fastq > barcode$i.fastq && pigz * && rm FA* && cd ..; done
