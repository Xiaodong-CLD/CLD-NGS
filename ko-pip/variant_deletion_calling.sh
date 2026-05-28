#!/bin/bash

# Set variables from command-line arguments
REFERENCE_PATH=$1
SORTED_BAM_PATH=$2
RESULT_FOLDER=$3


# Run Delly
delly call -g $REFERENCE_PATH -o $RESULT_FOLDER/delly.bcf $SORTED_BAM_PATH

# Run BCFtools for SNV calling and indexing
bcftools mpileup -Ou -f $REFERENCE_PATH $SORTED_BAM_PATH | bcftools call -mv -Oz -o $RESULT_FOLDER/variants.vcf.gz
bcftools index $RESULT_FOLDER/variants.vcf.gz

echo "[Takeda Pipeline Message] Processing completed. Check the results in $RESULT_FOLDER."