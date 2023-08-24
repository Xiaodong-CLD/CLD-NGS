#!/bin/bash

# Command-line arguments
main_dir=$1
plasmid_genome_dir=$2
virus_genome_dir=$3

# Loop over all BAM files in the directory
for bam_path in ${main_dir}/data/*.bam
do
    # Extract the sample name from the BAM file path
    sample_name=$(basename ${bam_path} .bam)

    # Create a directory for the sample
    mkdir -p ${main_dir}/${sample_name}
    cd ${main_dir}/${sample_name}

    # Check if unaligned fastq files already exist
    if [[ ! -f "unaligned_${sample_name}_1.fastq" || ! -f "unaligned_${sample_name}_2.fastq" ]]; then
        # Convert the BAM file to unaligned fastq files
        # Note: This requires samtools and bedtools to be installed
        samtools fastq -1 unaligned_${sample_name}_1.fastq -2 unaligned_${sample_name}_2.fastq ${bam_path}
    else
        echo "Unaligned fastq files for ${sample_name} already exist, skipping conversion."
    fi

    # Align the unaligned fastq files against the plasmid database
    STAR --genomeDir ${plasmid_genome_dir} --readFilesIn unaligned_${sample_name}_1.fastq unaligned_${sample_name}_2.fastq --runThreadN 24 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR_plasmid_${sample_name}_
    
    samtools coverage STAR_plasmid_${sample_name}_Aligned.sortedByCoord.out.bam > ${sample_name}-plasmid-coverage.txt
    python /mnt/efs/home/zds0975/pipeline_test/voof/unmap_analysis/scripts/analysis/coverage.py ${sample_name}-plasmid-coverage.txt -o ${sample_name}-plasmid-sorted-coverage.txt

    # Extract unmapped reads
    samtools view -b -f 4 STAR_plasmid_${sample_name}_Aligned.sortedByCoord.out.bam > STAR_plasmid_${sample_name}_unmapped.bam
    samtools fastq -1 unaligned_to_plasmid_${sample_name}_1.fastq -2 unaligned_to_plasmid_${sample_name}_2.fastq STAR_plasmid_${sample_name}_unmapped.bam

    # Create a directory for the virus
    mkdir -p virus
    cd virus

    # Check if coverage.txt already exists
    if [[ -f "${sample_name}-coverage.txt" ]]; then
        echo "${sample_name}-coverage.txt already exists, skipping steps for ${sample_name} virus."
        cd ${main_dir}
        continue
    fi

    # Align the unaligned fastq files against the virus database
    STAR --genomeDir ${virus_genome_dir} --readFilesIn ../unaligned_to_plasmid_${sample_name}_1.fastq ../unaligned_to_plasmid_${sample_name}_2.fastq --runThreadN 24 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR_virus_${sample_name}_

    # Calculate coverage and sequencing depth for each viral sequence in the database
    samtools coverage STAR_virus_${sample_name}_Aligned.sortedByCoord.out.bam > ${sample_name}-virus-coverage.txt
    python /mnt/efs/home/zds0975/pipeline_test/voof/unmap_analysis/scripts/analysis/coverage.py ${sample_name}-virus-coverage.txt -o ${sample_name}-virus-sorted-coverage.txt

    # Go back to the main directory for the next iteration
    cd ${main_dir}
done
