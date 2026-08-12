#!/bin/bash

# Command-line arguments
bam_dir=$1
plasmid_genome_dir=$2
virus_genome_dir=$3
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Loop over all BAM files in the directory
for bam_path in ${bam_dir}/*.bam
do
    # Extract the sample name from the BAM file path
    sample_name=$(basename ${bam_path} .bam)

    # Create a directory for the sample
    mkdir -p ${bam_dir}/../'Result'/${sample_name}
    cd ${bam_dir}/../'Result'/${sample_name}

    # Extract unmapped reads to an "unmapped.bam" file
    if [[ ! -f "unmapped_${sample_name}.bam" ]]; then
        samtools view -b -f 4 ${bam_path} > unmapped_${sample_name}.bam
    else
        echo "Unmapped bam file for ${sample_name} already exists, skipping extraction."
    fi
    unmapped_reads=$(samtools view -c unmapped_${sample_name}.bam)
    echo "Unmapped reads: $unmapped_reads" >> read_counts.txt


    # Convert the "unmapped.bam" file to unaligned fastq files for paired and single-end reads
    if [[ ! -f "unaligned_${sample_name}_1.fastq" ]] || [[ ! -f "unaligned_${sample_name}_2.fastq" ]] || [[ ! -f "unaligned_${sample_name}_single.fastq" ]]; then
        samtools fastq -1 "unaligned_${sample_name}_1.fastq" -2 "unaligned_${sample_name}_2.fastq" -s "unaligned_${sample_name}_single.fastq" unmapped_${sample_name}.bam
    else
        echo "Unaligned fastq files for ${sample_name} already exist, skipping conversion."
    fi

    # Align paired-end reads against the plasmid database
    if [[ -f "unaligned_${sample_name}_1.fastq" ]] && [[ -f "unaligned_${sample_name}_2.fastq" ]]; then
        STAR --genomeDir "${plasmid_genome_dir}" --readFilesIn "unaligned_${sample_name}_1.fastq" "unaligned_${sample_name}_2.fastq" \
             --runThreadN 24 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
             --outFileNamePrefix "STAR_plasmid_${sample_name}_paired_"
    fi

    # Align single-end reads against the plasmid database
    if [[ -f "unaligned_${sample_name}_single.fastq" ]]; then
        STAR --genomeDir "${plasmid_genome_dir}" --readFilesIn "unaligned_${sample_name}_single.fastq" \
             --runThreadN 24 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
             --outFileNamePrefix "STAR_plasmid_${sample_name}_single_"
    fi
    # Generate coverage reports for plasmid alignments
    samtools coverage STAR_plasmid_${sample_name}_paired_Aligned.sortedByCoord.out.bam > "${sample_name}-plasmid-paired-coverage.txt"
    samtools coverage STAR_plasmid_${sample_name}_single_Aligned.sortedByCoord.out.bam > "${sample_name}-plasmid-single-coverage.txt"
    python "$SCRIPT_DIR/coverage.py" "${sample_name}-plasmid-paired-coverage.txt" -o "${sample_name}-plasmid-paired-sorted-coverage.txt"
    python "$SCRIPT_DIR/coverage.py" "${sample_name}-plasmid-single-coverage.txt" -o "${sample_name}-plasmid-single-sorted-coverage.txt"
    
    # Count reads mapped and unmapped to plasmids
    mapped_plasmid_reads_paired=$(samtools view -c -F 4 STAR_plasmid_${sample_name}_paired_Aligned.sortedByCoord.out.bam)
    unmapped_plasmid_reads_paired=$(samtools view -c -f 4 STAR_plasmid_${sample_name}_paired_Aligned.sortedByCoord.out.bam)
    mapped_plasmid_reads_single=$(samtools view -c -F 4 STAR_plasmid_${sample_name}_single_Aligned.sortedByCoord.out.bam)
    unmapped_plasmid_reads_single=$(samtools view -c -f 4 STAR_plasmid_${sample_name}_single_Aligned.sortedByCoord.out.bam)
    total_mapped_plasmid_reads=$((mapped_plasmid_reads_paired + mapped_plasmid_reads_single))
    total_unmapped_plasmid_reads=$((unmapped_plasmid_reads_paired + unmapped_plasmid_reads_single))


    echo "Mapped to plasmids (paired-end): $mapped_plasmid_reads_paired" >> read_counts.txt
    echo "Unmapped to plasmids (paired-end): $unmapped_plasmid_reads_paired" >> read_counts.txt
    echo "Mapped to plasmids (single-end): $mapped_plasmid_reads_single" >> read_counts.txt
    echo "Unmapped to plasmids (single-end): $unmapped_plasmid_reads_single" >> read_counts.txt
    echo "Total Mapped to Plasmids: $total_mapped_plasmid_reads" >> read_counts.txt
    echo "Total Unmapped to Plasmids: $total_unmapped_plasmid_reads" >> read_counts.txt

    # Extract unmapped reads from plasmid alignments
    samtools view -b -f 4 STAR_plasmid_${sample_name}_paired_Aligned.sortedByCoord.out.bam > STAR_plasmid_${sample_name}_paired_unmapped.bam
    samtools view -b -f 4 STAR_plasmid_${sample_name}_single_Aligned.sortedByCoord.out.bam > STAR_plasmid_${sample_name}_single_unmapped.bam

    # Convert the plasmid unmapped bam files to fastq for paired and single-end reads
    samtools fastq -1 "unaligned_to_plasmid_${sample_name}_1.fastq" -2 "unaligned_to_plasmid_${sample_name}_2.fastq" STAR_plasmid_${sample_name}_paired_unmapped.bam
    samtools fastq "unaligned_to_plasmid_${sample_name}_single.fastq" STAR_plasmid_${sample_name}_single_unmapped.bam

    # Align reads against the virus database
    if [[ -f "unaligned_to_plasmid_${sample_name}_1.fastq" ]] && [[ -f "unaligned_to_plasmid_${sample_name}_2.fastq" ]]; then
        STAR --genomeDir "${virus_genome_dir}" --readFilesIn "unaligned_to_plasmid_${sample_name}_1.fastq" "unaligned_to_plasmid_${sample_name}_2.fastq" \
             --runThreadN 24 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
             --outFileNamePrefix "STAR_virus_${sample_name}_paired_"
    fi

    if [[ -f "unaligned_to_plasmid_${sample_name}_single.fastq" ]]; then
        STAR --genomeDir "${virus_genome_dir}" --readFilesIn "unaligned_to_plasmid_${sample_name}_single.fastq" \
             --runThreadN 24 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
             --outFileNamePrefix "STAR_virus_${sample_name}_single_"
    fi
    
    # Generate coverage reports for virus alignments
    samtools coverage STAR_virus_${sample_name}_paired_Aligned.sortedByCoord.out.bam > "${sample_name}-virus-paired-coverage.txt"
    samtools coverage STAR_virus_${sample_name}_single_Aligned.sortedByCoord.out.bam > "${sample_name}-virus-single-coverage.txt"
    python "$SCRIPT_DIR/coverage.py" "${sample_name}-virus-paired-coverage.txt" -o "${sample_name}-virus-paired-sorted-coverage.txt"
    python "$SCRIPT_DIR/coverage.py" "${sample_name}-virus-single-coverage.txt" -o "${sample_name}-virus-single-sorted-coverage.txt"

    # Count reads mapped and unmapped to virus
    mapped_virus_reads_paired=$(samtools view -c -F 4 STAR_virus_${sample_name}_paired_Aligned.sortedByCoord.out.bam)
    unmapped_virus_reads_paired=$(samtools view -c -f 4 STAR_virus_${sample_name}_paired_Aligned.sortedByCoord.out.bam)
    mapped_virus_reads_single=$(samtools view -c -F 4 STAR_virus_${sample_name}_single_Aligned.sortedByCoord.out.bam)
    unmapped_virus_reads_single=$(samtools view -c -f 4 STAR_virus_${sample_name}_single_Aligned.sortedByCoord.out.bam)
    total_mapped_virus_reads=$((mapped_virus_reads_paired + mapped_virus_reads_single))
    total_unmapped_virus_reads=$((unmapped_virus_reads_paired + unmapped_virus_reads_single))


    echo "Mapped to virus (paired-end): $mapped_virus_reads_paired" >> read_counts.txt
    echo "Unmapped to virus (paired-end): $unmapped_virus_reads_paired" >> read_counts.txt
    echo "Mapped to virus (single-end): $mapped_virus_reads_single" >> read_counts.txt
    echo "Unmapped to virus (single-end): $unmapped_virus_reads_single" >> read_counts.txt
    echo "Total Mapped to Virus: $total_mapped_virus_reads" >> read_counts.txt
    echo "Total Unmapped to Virus: $total_unmapped_virus_reads" >> read_counts.txt

    # Go back to the main directory for the next iteration
    cd ${bam_dir}
done
