#!/bin/bash

# Set variables from command-line arguments
REFERENCE_PATH=$1
SORTED_BAM_PATH=$2
OUTPUT_DIR=$3


# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Adding read group using Picard
OUTPUT_WITH_RG="$OUTPUT_DIR/output_with_RG.bam"
picard AddOrReplaceReadGroups \
      -I $SORTED_BAM_PATH \
      -O $OUTPUT_WITH_RG \
      -RGID 4 \
      -RGLB "lib1" \
      -RGPL "illumina" \
      -RGPU "unit1" \
      -RGSM "20" # read group sample name

SORTED_BAM_WITH_RG="$OUTPUT_DIR/sorted_output_with_RG.bam"
samtools sort -@ 16 $OUTPUT_WITH_RG -o $SORTED_BAM_WITH_RG
samtools index -@ 16 $SORTED_BAM_WITH_RG

# Update the SORTED_BAM_PATH variable to point to the BAM with added read group
SORTED_BAM_PATH=$SORTED_BAM_WITH_RG

# Extract the number of mapped reads from the BAM file
NUMBER_OF_READS=$(samtools view -c -F 4 $SORTED_BAM_PATH)
echo "[Takeda Pipeline] Number of mapped reads: $NUMBER_OF_READS"

# Create reference.dict file
gatk CreateSequenceDictionary -R $REFERENCE_PATH

Reads_Consider=$((NUMBER_OF_READS*2))

# Run GATK HaplotypeCaller for short variant discovery
HAPLOTYPE_OUTPUT_VCF="$OUTPUT_DIR/gatk_variants.vcf.gz"
gatk HaplotypeCaller \
    -R $REFERENCE_PATH \
    -I $SORTED_BAM_PATH \
    --max-reads-per-alignment-start $Reads_Consider \
    --stand-call-conf 20 \
    --active-probability-threshold 0.0002 \
    --min-base-quality-score 5 \
    -O $HAPLOTYPE_OUTPUT_VCF


if [ $? -eq 0 ]; then
    echo "[GATK Message] HaplotypeCaller executed successfully."
else
    echo "[GATK Message] HaplotypeCaller execution failed."
    exit 1
fi

echo "[GATK Message] Processing completed. Check the results in $OUTPUT_DIR."