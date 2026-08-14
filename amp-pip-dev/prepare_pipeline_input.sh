#!/bin/bash
# prepare_pipeline_input.sh -- Populate an input directory for the amp-pip-dev pipeline.
#
# Usage:
#   prepare_pipeline_input.sh <fastq_dir> <sample_sheet> <reference_file> <target_input_dir>
#
# Expected layout created under <target_input_dir>:
#   Fastq/         -- absolute symlinks to every FASTQ in <fastq_dir> (top-level)
#   sample_xlsx/   -- absolute symlink to <sample_sheet>
#   reference/     -- absolute symlink to <reference_file>
#
# Supported FASTQ extensions: .fastq .fastq.gz .fq .fq.gz
# Safe to re-run (ln -sfn); fails if a regular (non-symlink) file would be
# overwritten. Does not copy, move, or delete source files.
set -euo pipefail

usage() {
    echo "Usage: $(basename "$0") <fastq_dir> <sample_sheet> <reference_file> <target_input_dir>" >&2
    exit 1
}

[[ $# -eq 4 ]] || usage

FASTQ_DIR="$(realpath "$1")"
SAMPLE_SHEET="$(realpath "$2")"
REF_FILE="$(realpath "$3")"
TARGET_DIR="$(realpath -m "$4")"

# ── Validate sources ──────────────────────────────────────────────────────────
[[ -d "$FASTQ_DIR" ]]    || { echo "Error: fastq_dir not found: $FASTQ_DIR" >&2; exit 1; }
[[ -f "$SAMPLE_SHEET" ]] || { echo "Error: sample_sheet not found: $SAMPLE_SHEET" >&2; exit 1; }
[[ -f "$REF_FILE" ]]     || { echo "Error: reference_file not found: $REF_FILE" >&2; exit 1; }

# ── Create target structure ───────────────────────────────────────────────────
FASTQ_TARGET="${TARGET_DIR}/Fastq"
SHEET_TARGET="${TARGET_DIR}/sample_xlsx"
REF_TARGET="${TARGET_DIR}/reference"
mkdir -p "$FASTQ_TARGET" "$SHEET_TARGET" "$REF_TARGET"

# ── Symlink FASTQ files (top-level only; pipeline reads filenames from sheet) ─
n_fq=0
for src in \
    "${FASTQ_DIR}"/*.fastq \
    "${FASTQ_DIR}"/*.fastq.gz \
    "${FASTQ_DIR}"/*.fq \
    "${FASTQ_DIR}"/*.fq.gz
do
    [[ -e "$src" ]] || continue
    dest="${FASTQ_TARGET}/$(basename "$src")"
    if [[ -e "$dest" && ! -L "$dest" ]]; then
        echo "Error: regular file already occupies $dest" >&2; exit 1
    fi
    ln -sfn "$src" "$dest"
    n_fq=$(( n_fq + 1 ))
done

[[ $n_fq -gt 0 ]] || {
    echo "Error: no FASTQ files (.fastq/.fastq.gz/.fq/.fq.gz) found in $FASTQ_DIR" >&2
    exit 1
}

# ── Symlink sample sheet into sample_xlsx/ ────────────────────────────────────
SHEET_DEST="${SHEET_TARGET}/$(basename "$SAMPLE_SHEET")"
if [[ -e "$SHEET_DEST" && ! -L "$SHEET_DEST" ]]; then
    echo "Error: regular file already occupies $SHEET_DEST" >&2; exit 1
fi
ln -sfn "$SAMPLE_SHEET" "$SHEET_DEST"

# ── Symlink reference file ───────────────────────────────────────────────────
REF_DEST="${REF_TARGET}/$(basename "$REF_FILE")"
if [[ -e "$REF_DEST" && ! -L "$REF_DEST" ]]; then
    echo "Error: regular file already occupies $REF_DEST" >&2; exit 1
fi
ln -sfn "$REF_FILE" "$REF_DEST"

# ── Summary ───────────────────────────────────────────────────────────────────
echo "Input directory populated: $TARGET_DIR"
echo "  sample sheet : $(basename "$SHEET_DEST") -> $SAMPLE_SHEET"
echo "  FASTQ links  : $n_fq file(s) in $FASTQ_TARGET"
echo "  ref file     : $(basename "$REF_DEST") -> $REF_FILE"
