#!/bin/bash

mkdir consensus
cd consensus
for filename in ../bc_split/*.fasta; do
    echo "Generating consensus reads for: $filename"
    python ../generate_consensus.py -i "$filename" -o "$(basename "$filename" .fasta)_consensus.fasta"
done
cd ..

mkdir error_analysis
cd error_analysis

for filename in ../consensus/*.fasta; do
    bwa index ../template.fasta
    echo "Calculating errors for: $filename"
    bwa aln -n 10 -o 10 -O 0 -E 0 -t 8 ../template.fasta "$filename" > aligned.sai
    bwa samse ../template.fasta aligned.sai "$filename" > aligned.sam
    samtools sort -T /tmp/aln.sorted -o aligned-sorted.bam aligned.sam
    samtools faidx ../template.fasta
    samtools mpileup -d 10000000 -f ../template.fasta aligned-sorted.bam > aligned-sorted.pileup
    python ../error_analysis_pileup.py -pileup aligned-sorted.pileup > "$(basename "$filename" .fasta)_err.tsv"

done

