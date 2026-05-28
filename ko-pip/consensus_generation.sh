#!/bin/bash

# Define arguments
reference_fa=$1
alignments_bam=$2
output_dir=$3
non_deleted_regions=$4

# Define directory for storing intermediate files
temp_dir="${output_dir}/temp"

# Create output directory if it doesn't exist
mkdir -p "${output_dir}"
mkdir -p "${temp_dir}"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Process each non-deleted region
while read -r start end; do
    # Define region_string
    region_string="${start}-${end}"
    
    # Extract region from reference genome
    python3 "$SCRIPT_DIR/extract_subsequence.py" "${reference_fa}" "${start}" "${end}" "${temp_dir}/${region_string}.fa"

    # Skip to next iteration if consensus sequence already exists
    if [ -f "${output_dir}/${region_string}_consensus.fa" ]; then
        echo "[Takeda Pipeline Message] Consensus sequence for region ${region_string} already exists. Skipping..."
        continue
    fi

    # Call variants
    bcftools mpileup -Ou -f "${temp_dir}/${region_string}.fa" "${alignments_bam}" | bcftools call -mv -Oz -o "${temp_dir}/${region_string}_calls.vcf.gz"
    bcftools index "${temp_dir}/${region_string}_calls.vcf.gz"

    # Normalize indels
    bcftools norm -f "${temp_dir}/${region_string}.fa" "${temp_dir}/${region_string}_calls.vcf.gz" -Ob -o "${temp_dir}/${region_string}_calls.norm.bcf"

    # Filter adjacent indels within 5bp
    bcftools filter --IndelGap 5 "${temp_dir}/${region_string}_calls.norm.bcf" -Ob -o "${temp_dir}/${region_string}_calls.norm.flt-indels.bcf"

    # Apply variants to create consensus sequence
    cat "${temp_dir}/${region_string}.fa" | bcftools consensus "${temp_dir}/${region_string}_calls.vcf.gz" > "${output_dir}/${region_string}_consensus.fa"
done < "${non_deleted_regions}"

# Concatenate all consensus sequences
cat ${output_dir}/*_consensus.fa > "${output_dir}/consensus.fa"
