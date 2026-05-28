#!/bin/bash

# Usage: bash gene_annotation_pipeline.sh [-i gene_of_interest.fasta] [-r human_genome_reference.fasta] [-a human_annotation.gtf] [-o output_directory]

# Default file paths and output directory
gene_of_interest=""
human_genome_reference=""
human_annotation=""
output_directory=""

# Function to display usage instructions
function display_usage {
    echo "Usage: bash gene_annotation_pipeline.sh [-i gene_of_interest.fasta] [-r human_genome_reference.fasta] [-a human_annotation.gtf] [-o output_directory]"
    echo "Options:"
    echo "  -i, --gene-of-interest   Path to the gene of interest sequence (FASTA)"
    echo "  -r, --genome-reference   Path to the human genome reference sequence (FASTA)"
    echo "  -a, --annotation         Path to the human annotation file (GTF)"
    echo "  -o, --output-directory   Directory to store the output files"
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--gene-of-interest)
            gene_of_interest="$2"
            shift 2 ;;
        -r|--genome-reference)
            human_genome_reference="$2"
            shift 2 ;;
        -a|--annotation)
            human_annotation="$2"
            shift 2 ;;
        -o|--output-directory)
            output_directory="$2"
            shift 2 ;;
        *)
            display_usage
            exit 1 ;;
    esac
done

# Check if any required input is missing
if [[ -z $gene_of_interest || -z $human_genome_reference || -z $human_annotation ]]; then
    echo "Error: Missing required input. Please provide the gene of interest sequence, human genome reference, and human annotation file."
    display_usage
    exit 1
fi

# Set the output directory
if [[ -z $output_directory ]]; then
    output_directory=$(pwd)
fi

# Rest of the pipeline steps

# Step 1: Align the gene of interest to the human genome reference using Minimap2
minimap2 -ax map-ont -t 8 ${human_genome_reference} ${gene_of_interest} > ${output_directory}/gene_of_interest.sam

# Step 2: Convert the SAM file to BAM format
samtools view -bS ${output_directory}/gene_of_interest.sam > ${output_directory}/gene_of_interest.bam

# Step 3: Sort the BAM file
samtools sort -o ${output_directory}/gene_of_interest.sorted.bam ${output_directory}/gene_of_interest.bam

# Step 4: Index the sorted BAM file
samtools index ${output_directory}/gene_of_interest.sorted.bam

# Step 5: Perform gene annotation using StringTie
stringtie -G ${human_annotation} -e -B -p 8 -o ${output_directory}/gene_of_interest_annotation.gtf ${output_directory}/gene_of_interest.sorted.bam

# Step 6: Optional - Merge transcripts from multiple samples (if applicable)
# If you have multiple samples, you can merge the individual StringTie output using the following command:
# stringtie --merge -G ${human_annotation} -o merged_transcripts.gtf merged_list.txt

# Step 7: Optional - Refine the merged transcript annotations (if applicable)
# If you performed the previous step, you can refine the merged transcript annotations using the following command:
# stringtie -eB -p 8 -G merged_transcripts.gtf -o refined_transcripts.gtf gene_of_interest.sorted.bam

echo "Gene annotation for the gene of interest has been performed."
