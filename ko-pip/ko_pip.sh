#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> -r <reference_fa> [-s <session_name>]"
    echo
    echo "  -i  Input directory (must contain sample_xlsx/ and Fastq/)"
    echo "  -o  Output directory"
    echo "  -r  Reference fasta"
    echo "  -s  tmux session name (default: derived from output dir basename)"
    exit 1
}

INPUT=""
OUTPUT=""
REFERENCE=""
SESSION=""

while getopts ":i:o:r:s:h" opt; do
    case "$opt" in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        s) SESSION="$OPTARG" ;;
        h|*) usage ;;
    esac
done

if [[ -z "$INPUT" || -z "$OUTPUT" || -z "$REFERENCE" ]]; then
    usage
fi

if [[ -z "$SESSION" ]]; then
    SESSION="$(basename "${OUTPUT%/}")"
fi

mkdir -p logs "$OUTPUT"
LOG="logs/${SESSION}.log"

echo "Starting tmux session '$SESSION'"
echo "  input:     $INPUT"
echo "  output:    $OUTPUT"
echo "  reference: $REFERENCE"
echo "  log:       $LOG"

tmux new-session -d -s "$SESSION" \
    "python3 -u pipeline_combined.py \
        -i '$INPUT' \
        -o '$OUTPUT' \
        -r '$REFERENCE' \
        > '$LOG' 2>&1"

echo "Tail the log with: tail -f $LOG"
echo "Attach with:       tmux attach -t $SESSION"
