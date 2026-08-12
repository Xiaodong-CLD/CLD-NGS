#!/bin/bash
set -euo pipefail

CONDA_ENV="unmap_analysis"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PIPELINE_SCRIPT="${SCRIPT_DIR}/scripts/analysis/unmap_analysis.sh"

usage() {
    echo "Usage: $0 -b <bam_dir> -p <plasmid_fa_dir> -v <virus_genome_dir> [-s <session_name>]"
    echo
    echo "  -b  Directory containing input .bam files"
    echo "  -p  Directory containing plasmid .fa files (STAR index will be built here)"
    echo "  -v  Pre-built STAR virus genome directory"
    echo "  -s  tmux session name (default: derived from bam_dir basename)"
    exit 1
}

BAM_DIR=""
PLASMID_FA_DIR=""
VIRUS_GENOME_DIR=""
SESSION=""

while getopts ":b:p:v:s:h" opt; do
    case "$opt" in
        b) BAM_DIR="$OPTARG" ;;
        p) PLASMID_FA_DIR="$OPTARG" ;;
        v) VIRUS_GENOME_DIR="$OPTARG" ;;
        s) SESSION="$OPTARG" ;;
        h|*) usage ;;
    esac
done

if [[ -z "$BAM_DIR" || -z "$PLASMID_FA_DIR" || -z "$VIRUS_GENOME_DIR" ]]; then
    usage
fi

if [[ -z "$SESSION" ]]; then
    SESSION="$(basename "${BAM_DIR%/}")"
fi

# ── 1. Dependency check ───────────────────────────────────────────────────────

CONDA_BASE="$(conda info --base 2>/dev/null)"
if [[ -z "$CONDA_BASE" || ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
    echo "Error: could not locate conda installation (is conda on PATH?)" >&2
    exit 1
fi

if ! conda env list | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
    echo "Error: conda env '$CONDA_ENV' not found. Create it with:" >&2
    echo "  conda env create -f envs/unmap_analysis.yml" >&2
    exit 1
fi

for tool in samtools STAR; do
    if ! conda run -n "$CONDA_ENV" which "$tool" &>/dev/null; then
        echo "Error: '$tool' not found in conda env '$CONDA_ENV'" >&2
        exit 1
    fi
done

# ── 2. Setup ──────────────────────────────────────────────────────────────────

mkdir -p "${BAM_DIR}/../Result" "${BAM_DIR}/../logs"
LOG="${BAM_DIR}/../logs/${SESSION}.log"
PLASMID_GENOME_DIR="${PLASMID_FA_DIR}/plasmid_star_dir"

echo "Starting tmux session '$SESSION'"
echo "  bam_dir:           $BAM_DIR"
echo "  plasmid_fa_dir:    $PLASMID_FA_DIR"
echo "  plasmid_genome_dir:$PLASMID_GENOME_DIR"
echo "  virus_genome_dir:  $VIRUS_GENOME_DIR"
echo "  env:               $CONDA_ENV"
echo "  log:               $LOG"

# ── 3. Run pipeline in tmux session ──────────────────────────────────────────
# Mirrors job.script logic: concatenate plasmid FAs → build STAR index → run pipeline

tmux new-session -d -s "$SESSION" \
    "source '$CONDA_BASE/etc/profile.d/conda.sh' && \
     conda activate '$CONDA_ENV' && \
     mkdir -p '$PLASMID_GENOME_DIR' && \
     cat '${PLASMID_FA_DIR}'/*.fa > '${PLASMID_FA_DIR}/plasmid.fa' && \
     STAR --runMode genomeGenerate \
          --genomeDir '$PLASMID_GENOME_DIR' \
          --genomeFastaFiles '${PLASMID_FA_DIR}/plasmid.fa' \
          --runThreadN 24 && \
     bash '$PIPELINE_SCRIPT' '$BAM_DIR' '$PLASMID_GENOME_DIR' '$VIRUS_GENOME_DIR' \
     > '$LOG' 2>&1"

echo "Tail the log with: tail -f $LOG"
echo "Attach with:       tmux attach -t $SESSION"
