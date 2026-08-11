#!/usr/bin/env bash
# rna_var.sh — local runner for the RNAseq-var variant-calling pipeline.
#
# Checks dependencies, organizes inputs, and runs pipeline_run.py inside a
# detached tmux session.
#
# Usage
# -----
#   bash rna_var.sh \
#       --fastq-dir   /path/to/Fastq \
#       --sample-sheet /path/to/samples.xlsx \
#       --reference   /path/to/reference.fna \
#       --goi         /path/to/goi.fa \
#       --run-name    my_run_name \
#       [--run-dir    /custom/output/root]   # default: ./runs/<run_name>
#
# The tmux session is named  rnaseq_var_<run_name>.
# Attach to it with:
#   tmux attach -t rnaseq_var_<run_name>
#
# Pipeline outputs are written to:
#   <run_dir>/Result/

set -euo pipefail

# ---------------------------------------------------------------------------
# Script location (used to find sibling Python scripts)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
FASTQ_DIR=""
SAMPLE_SHEET=""
REFERENCE=""
GOI=""
RUN_NAME=""
RUN_DIR=""

usage() {
    sed -n '3,20p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq-dir)    FASTQ_DIR="$(realpath "$2")";    shift 2 ;;
        --sample-sheet) SAMPLE_SHEET="$(realpath "$2")"; shift 2 ;;
        --reference)    REFERENCE="$(realpath "$2")";    shift 2 ;;
        --goi)          GOI="$(realpath "$2")";          shift 2 ;;
        --run-name)     RUN_NAME="$2";                   shift 2 ;;
        --run-dir)      RUN_DIR="$(realpath "$2")";      shift 2 ;;
        -h|--help)      usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# Required argument checks
missing=0
for flag_name in FASTQ_DIR SAMPLE_SHEET REFERENCE GOI RUN_NAME; do
    if [[ -z "${!flag_name}" ]]; then
        echo "[rna_var] ERROR: --$(echo "$flag_name" | tr '[:upper:]_' '[:lower:]-') is required."
        missing=1
    fi
done
[[ $missing -eq 1 ]] && usage

# Default run directory
if [[ -z "$RUN_DIR" ]]; then
    RUN_DIR="$SCRIPT_DIR/runs/$RUN_NAME"
fi

# ---------------------------------------------------------------------------
# Dependency checks
# ---------------------------------------------------------------------------
echo "[rna_var] Checking dependencies..."

check_cmd() {
    if ! command -v "$1" &>/dev/null; then
        echo "[rna_var] ERROR: '$1' not found in PATH. Activate the RNAseq conda environment first."
        return 1
    fi
    echo "[rna_var]   OK  $1"
}

DEPS_OK=1
for tool in python conda tmux bwa samtools fastp varscan pindel Rscript; do
    check_cmd "$tool" || DEPS_OK=0
done

# Check Python packages
python -c "import pandas, openpyxl" 2>/dev/null \
    && echo "[rna_var]   OK  python: pandas, openpyxl" \
    || { echo "[rna_var] ERROR: Python packages 'pandas' or 'openpyxl' missing."; DEPS_OK=0; }

if [[ $DEPS_OK -eq 0 ]]; then
    echo ""
    echo "[rna_var] Dependency check failed. Run the following and retry:"
    echo "    conda activate RNAseq"
    exit 1
fi
echo "[rna_var] All dependencies satisfied."
echo ""

# ---------------------------------------------------------------------------
# Validate input paths
# ---------------------------------------------------------------------------
echo "[rna_var] Validating input paths..."
error=0

[[ -d "$FASTQ_DIR" ]]     || { echo "[rna_var] ERROR: FASTQ directory not found: $FASTQ_DIR";     error=1; }
[[ -f "$SAMPLE_SHEET" ]]  || { echo "[rna_var] ERROR: Sample sheet not found: $SAMPLE_SHEET";     error=1; }
[[ -f "$REFERENCE" ]]     || { echo "[rna_var] ERROR: Reference FASTA not found: $REFERENCE";     error=1; }
[[ -f "$GOI" ]]           || { echo "[rna_var] ERROR: GOI FASTA not found: $GOI";                 error=1; }
[[ $error -eq 0 ]]        || exit 1

echo "[rna_var] Inputs OK."
echo ""

# ---------------------------------------------------------------------------
# Organize input directory structure (symlinks)
# ---------------------------------------------------------------------------
echo "[rna_var] Organizing input structure in: $RUN_DIR"
mkdir -p "$RUN_DIR"

python "$SCRIPT_DIR/organize_inputs.py" \
    --fastq-dir    "$FASTQ_DIR" \
    --sample-sheet "$SAMPLE_SHEET" \
    --reference    "$REFERENCE" \
    --goi          "$GOI" \
    --run-dir      "$RUN_DIR" \
    --overwrite
echo ""

# ---------------------------------------------------------------------------
# Log file
# ---------------------------------------------------------------------------
LOG_FILE="$RUN_DIR/pipeline.log"
touch "$LOG_FILE"
echo "[rna_var] Pipeline log: $LOG_FILE"

# ---------------------------------------------------------------------------
# Launch pipeline in a detached tmux session
# ---------------------------------------------------------------------------
SESSION_NAME="rnaseq_var_${RUN_NAME}"

# Kill any existing session with the same name so re-runs are clean.
tmux kill-session -t "$SESSION_NAME" 2>/dev/null || true

echo "[rna_var] Starting tmux session: $SESSION_NAME"
echo "[rna_var] Attach with:  tmux attach -t $SESSION_NAME"
echo ""

# The pipeline reads sample_xlsx/*.xlsx relative to CWD, so cd into RUN_DIR.
tmux new-session -d -s "$SESSION_NAME" \
    "cd '$RUN_DIR' && python -u '$SCRIPT_DIR/pipeline_run.py' '$RUN_DIR' 2>&1 | tee '$LOG_FILE'; \
     echo '[rna_var] Pipeline finished with exit code '\$?'. Press any key to close.'; read -n1"

echo "[rna_var] Pipeline launched. Monitor progress:"
echo "    tail -f $LOG_FILE"
echo "    tmux attach -t $SESSION_NAME"
