#!/usr/bin/env python3
"""User-facing pipeline launcher for amp-pip-dev.

Accepts the four user-visible inputs and runs the full amplicon
variant-calling pipeline:

    python run_pipeline.py \
        --sample-sheet  samples.xlsx \
        --fastq-dir     /path/to/Fastq \
        --reference     goi.fa \
        --output-dir    /path/to/output

How it works
------------
1. tools/prepare_pipeline_input.sh creates the directory structure that
   pipeline_combined.py expects inside <output-dir>:

       <output-dir>/Fastq/           <- symlinks to source FASTQ files
       <output-dir>/sample_xlsx/     <- symlink to sample sheet
       <output-dir>/reference/       <- symlink to reference FASTA

2. pipeline_combined.py is then called with <output-dir> as its input_dir.
   It writes Results/ and tmp/ there.

Note: this is a compatibility shim.  See docs/upgrade-2026-08.md for the
planned refactoring that will allow pipeline_combined.py to accept these
paths directly and remove the symlink preparation step.
"""

import argparse
import pathlib
import subprocess
import sys

_HERE = pathlib.Path(__file__).parent.resolve()
_PREPARE_SH = _HERE / "tools" / "prepare_pipeline_input.sh"
_PIPELINE_COMBINED = _HERE / "pipeline_combined.py"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the amp-pip-dev amplicon variant-calling pipeline.",
    )
    parser.add_argument(
        "--sample-sheet", required=True,
        help="Path to the sample sheet (.xlsx).",
    )
    parser.add_argument(
        "--fastq-dir", required=True,
        help="Directory containing FASTQ files.",
    )
    parser.add_argument(
        "--reference", required=True,
        help="Amplicon reference FASTA file.",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Output directory (created if it does not exist).",
    )
    args = parser.parse_args()

    output_dir = pathlib.Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: populate <output-dir> with symlinks to user-supplied inputs.
    subprocess.run(
        [
            "bash", str(_PREPARE_SH),
            args.fastq_dir,
            args.sample_sheet,
            args.reference,
            str(output_dir),
        ],
        check=True,
    )

    # Step 2: run the pipeline on the prepared directory.
    subprocess.run(
        [sys.executable, str(_PIPELINE_COMBINED), str(output_dir)],
        check=True,
    )


if __name__ == "__main__":
    main()
