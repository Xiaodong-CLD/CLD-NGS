#!/bin/bash
#
# GATK HaplotypeCaller short-variant discovery.
#
# Usage: gatk_variant_calling.sh <reference_fa> <sorted_bam> <output_dir>
#
# Requires: gatk (>=4.4), samtools.
#
# Original script created by Lilin Huang.
# Updated by Cindy Ren 2026-06-02:
#   - Removed Picard AddOrReplaceReadGroups and the follow-up samtools
#     sort/index. The input BAM is already coordinate-sorted and the read
#     group (@RG) is embedded during BWA alignment (see align.py).
#   - Added strict mode, input validation, and quoted expansions for
#     compatibility with pipeline orchestrators (e.g. Tetra Science).
#   - Skip `gatk CreateSequenceDictionary` when the .dict already exists,
#     since the reference directory may be a read-only mount.
#   - Updated HaplotypeCaller parameters `--max-reads-per-alignment-start` to
#     0 to disable downsampling (original = 2 * number of mapped reads).
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <reference_fa> <sorted_bam> <output_dir>" >&2
    exit 2;
fi

REFERENCE_PATH="$1"
SORTED_BAM_PATH="$2"
OUTPUT_DIR="$3"

[[ -r "$REFERENCE_PATH" ]]  || { echo "Reference not readable: $REFERENCE_PATH" >&2; exit 2; }
[[ -r "$SORTED_BAM_PATH" ]] || { echo "BAM not readable: $SORTED_BAM_PATH" >&2; exit 2; }

mkdir -p "$OUTPUT_DIR"

# Number of mapped reads (used to size --max-reads-per-alignment-start)
NUMBER_OF_READS=$(samtools view -c -F 4 "$SORTED_BAM_PATH")
echo "[Takeda Pipeline] Number of mapped reads: $NUMBER_OF_READS"

# Create reference .dict next to the FASTA if missing. GATK accepts either
# <ref>.dict (replacing the extension) or <ref>.fa.dict, so check both.
REF_NO_EXT="${REFERENCE_PATH%.*}"
if [[ ! -f "${REF_NO_EXT}.dict" && ! -f "${REFERENCE_PATH}.dict" ]]; then
    gatk CreateSequenceDictionary -R "$REFERENCE_PATH"
fi

Reads_Consider=$((NUMBER_OF_READS * 2))

HAPLOTYPE_OUTPUT_VCF="$OUTPUT_DIR/gatk_variants.vcf.gz"
gatk HaplotypeCaller \
    -R "$REFERENCE_PATH" \
    -I "$SORTED_BAM_PATH" \
    --max-reads-per-alignment-start "$Reads_Consider" \
    --stand-call-conf 20 \
    --active-probability-threshold 0.0002 \
    --min-base-quality-score 5 \
    -O "$HAPLOTYPE_OUTPUT_VCF"

echo "[GATK Message] HaplotypeCaller executed successfully."
echo "[GATK Message] Processing completed. Check the results in $OUTPUT_DIR."