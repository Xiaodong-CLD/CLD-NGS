#!/bin/bash
set -euo pipefail

CONDA_ENV="ko_pip"

usage() {
    echo "Usage: $0 -i <fastq_dir> -p <sample_sheet> -o <output_dir> -r <reference_fa> -b <blast_ref_dir> [-s <session_name>]"
    echo
    echo "  -i  Directory containing fastq files"
    echo "  -p  Path to the sample sheet (.xlsx)"
    echo "  -o  Output directory"
    echo "  -r  Reference fasta"
    echo "  -b  Directory of plasmid FASTA files (*.fa/*.fasta) to BLAST breakpoint clusters against (required)"
    echo "  -s  tmux session name (default: derived from output dir basename)"
    exit 1
}

INPUT=""
SAMPLE_SHEET=""
OUTPUT=""
REFERENCE=""
SESSION=""
BLAST_REF=""

while getopts ":i:p:o:r:b:s:h" opt; do
    case "$opt" in
        i) INPUT="$OPTARG" ;;
        p) SAMPLE_SHEET="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        b) BLAST_REF="$OPTARG" ;;
        s) SESSION="$OPTARG" ;;
        h|*) usage ;;
    esac
done

if [[ -z "$INPUT" || -z "$SAMPLE_SHEET" || -z "$OUTPUT" || -z "$REFERENCE" || -z "$BLAST_REF" ]]; then
    usage
fi

if [[ -z "$SESSION" ]]; then
    SESSION="$(basename "${OUTPUT%/}")"
fi

mkdir -p logs "$OUTPUT"
LOG="logs/${SESSION}.log"

# Locate conda so `conda activate` works inside the non-interactive tmux shell
CONDA_BASE="$(conda info --base 2>/dev/null)"
if [[ -z "$CONDA_BASE" || ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
    echo "Error: could not locate conda installation (is conda on PATH?)" >&2
    exit 1
fi

if ! conda env list | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
    echo "Error: conda env '$CONDA_ENV' not found. Create it first." >&2
    exit 1
fi

echo "Starting tmux session '$SESSION'"
echo "  fastq_dir:    $INPUT"
echo "  sample_sheet: $SAMPLE_SHEET"
echo "  output:       $OUTPUT"
echo "  reference:    $REFERENCE"
echo "  blast_ref: $BLAST_REF"
echo "  env:       $CONDA_ENV"
echo "  log:       $LOG"

tmux new-session -d -s "$SESSION" \
    "source '$CONDA_BASE/etc/profile.d/conda.sh' && \
     conda activate '$CONDA_ENV' && \
     python3 -u pipeline.py \
        -i '$INPUT' \
        -s '$SAMPLE_SHEET' \
        -o '$OUTPUT' \
        -r '$REFERENCE' \
        --blast_ref '$BLAST_REF' \
        > '$LOG' 2>&1"

echo "Tail the log with: tail -f $LOG"
echo "Attach with:       tmux attach -t $SESSION"
