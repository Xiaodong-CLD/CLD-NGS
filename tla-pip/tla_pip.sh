#!/bin/bash
set -euo pipefail

CONDA_ENV="tla_pip"

usage() {
    echo "Usage: $0 -s <sample_sheet> -f <fastq_dir> -v <vector_fasta> -o <output_dir> [-r <references_dir>] [-a <accession>] [-n <session_name>] [--force-download]"
    echo
    echo "  -s  Sample sheet (.xlsx or .tsv)"
    echo "  -f  Directory containing FASTQ files"
    echo "  -v  Vector/GOI FASTA"
    echo "  -o  Output directory"
    echo "  -r  Local references directory (host genome + mixed ref cache; dev/testing override)"
    echo "  -a  Host genome NCBI accession (default: GCF_003668045.3)"
    echo "  -n  tmux session name (default: derived from output dir basename)"
    echo "  --force-download  Re-download host genome even if a cached copy exists"
    exit 1
}

SAMPLE_SHEET=""
FASTQ_DIR=""
VECTOR_FASTA=""
OUTPUT=""
REFERENCES=""
ACCESSION="GCF_003668045.3"
SESSION=""
FORCE_DOWNLOAD=""

while getopts ":s:f:v:o:r:a:n:h-:" opt; do
    case "$opt" in
        s) SAMPLE_SHEET="$OPTARG" ;;
        f) FASTQ_DIR="$OPTARG" ;;
        v) VECTOR_FASTA="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        r) REFERENCES="$OPTARG" ;;
        a) ACCESSION="$OPTARG" ;;
        n) SESSION="$OPTARG" ;;
        -) [[ "$OPTARG" == "force-download" ]] && FORCE_DOWNLOAD="--force-download" ;;
        h|*) usage ;;
    esac
done

if [[ -z "$SAMPLE_SHEET" || -z "$FASTQ_DIR" || -z "$VECTOR_FASTA" || -z "$OUTPUT" ]]; then
    usage
fi

if [[ -z "$SESSION" ]]; then
    SESSION="$(basename "${OUTPUT%/}")"
fi

mkdir -p "$OUTPUT/logs"
LOG="$OUTPUT/logs/${SESSION}.log"

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
echo "  sample sheet: $SAMPLE_SHEET"
echo "  fastq dir:    $FASTQ_DIR"
echo "  vector fasta: $VECTOR_FASTA"
echo "  output:       $OUTPUT"
echo "  references:   ${REFERENCES:-'(auto: <output>/reference/)'}"
echo "  accession:    $ACCESSION"
echo "  env:          $CONDA_ENV"
echo "  log:          $LOG"

# Build optional flags only if set
REF_ARG=""
[[ -n "$REFERENCES" ]] && REF_ARG="-r '$REFERENCES'"

tmux new-session -d -s "$SESSION" \
    "source '$CONDA_BASE/etc/profile.d/conda.sh' && \
     conda activate '$CONDA_ENV' && \
     python3 -u pipeline.py \
        -s '$SAMPLE_SHEET' \
        -f '$FASTQ_DIR' \
        -v '$VECTOR_FASTA' \
        -o '$OUTPUT' \
        -a '$ACCESSION' \
        $REF_ARG \
        $FORCE_DOWNLOAD \
        > '$LOG' 2>&1"

echo "Tail the log with: tail -f $LOG"
echo "Attach with:       tmux attach -t $SESSION"
