#!/bin/bash

# Check if the input file argument is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <input_bam>"
  exit 1
fi

# Input file
input_bam="$1"
input_directory=$(dirname "$input_bam")
input_filename=$(basename "${input_bam%.*}")

# Generate output file names
output_unmapped_bam="$input_directory/${input_filename}.unmapped.bam"
output_R1_fq="$input_directory/unmap_R1.fq"
output_R2_fq="$input_directory/unmap_R2.fq"
output_assembled_fasta="$input_directory/sorted.contigs.fasta"

# Step 1: Extract unmapped reads
samtools view -b -f 4 "$input_bam" > "$output_unmapped_bam"

# Step 2: Convert unmapped reads to FASTQ format
bedtools bamtofastq -i "$output_unmapped_bam" -fq "$output_R1_fq" -fq2 "$output_R2_fq"

# Step 3: Assemble contigs using SPAdes
spades.py -1 "$output_R1_fq" -2 "$output_R2_fq" --only-assembler -o "$input_directory" --sc --careful
