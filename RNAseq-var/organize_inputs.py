"""organize_inputs.py — Assemble the RNAseq-var input directory structure.

RNAseq-var/pipeline_run.py expects input files to be pre-arranged in four
subdirectories.  This script creates that layout inside a run directory using
symbolic links, so the originals are never copied or moved.

Expected output layout
----------------------
    <run_dir>/
        Fastq/
            <all fastq files from --fastq-dir>  (symlinks)
        reference/
            <reference FASTA filename>           (symlink)
        goi/
            <GOI FASTA filename>                 (symlink)
        sample_xlsx/
            <sample sheet filename>              (symlink)

Usage
-----
    python organize_inputs.py \\
        --fastq-dir    /path/to/raw_fastqs \\
        --sample-sheet /path/to/samples.xlsx \\
        --reference    /path/to/reference.fna \\
        --goi          /path/to/goi.fa \\
        --run-dir      /path/to/run_output_dir \\
        [--overwrite]
"""

from __future__ import annotations

import argparse
import pathlib
import sys


_FASTQ_EXTS = {".fastq", ".fq", ".gz"}


def _is_fastq(p: pathlib.Path) -> bool:
    return p.name.endswith(".fastq.gz") or p.name.endswith(".fq.gz") \
        or p.name.endswith(".fastq") or p.name.endswith(".fq")


def _make_symlink(src: pathlib.Path, dst: pathlib.Path, overwrite: bool) -> str:
    """Create dst → src symlink.  Returns 'created', 'skipped', or 'replaced'."""
    if dst.exists() or dst.is_symlink():
        if not overwrite:
            return "skipped"
        dst.unlink()
        dst.symlink_to(src)
        return "replaced"
    dst.symlink_to(src)
    return "created"


def organize(
    fastq_dir: pathlib.Path,
    sample_sheet: pathlib.Path,
    reference: pathlib.Path,
    goi: pathlib.Path,
    run_dir: pathlib.Path,
    overwrite: bool = False,
) -> None:
    """Create the four-subdirectory layout inside *run_dir*.

    Raises ValueError on missing inputs.
    Raises FileExistsError if a symlink already exists and overwrite=False.
    """
    if not fastq_dir.is_dir():
        raise ValueError(f"FASTQ directory not found: {fastq_dir}")
    if not sample_sheet.is_file():
        raise ValueError(f"Sample sheet not found: {sample_sheet}")
    if not reference.is_file():
        raise ValueError(f"Reference FASTA not found: {reference}")
    if not goi.is_file():
        raise ValueError(f"GOI FASTA not found: {goi}")

    fastq_files = sorted(f for f in fastq_dir.iterdir() if f.is_file() and _is_fastq(f))
    if not fastq_files:
        raise ValueError(f"No FASTQ files found in {fastq_dir}")

    # Create subdirectory structure.
    fastq_subdir = run_dir / "Fastq"
    ref_subdir   = run_dir / "reference"
    goi_subdir   = run_dir / "goi"
    xlsx_subdir  = run_dir / "sample_xlsx"
    for d in (fastq_subdir, ref_subdir, goi_subdir, xlsx_subdir):
        d.mkdir(parents=True, exist_ok=True)

    # Symlink FASTQ files.
    print(f"[organize_inputs] Linking {len(fastq_files)} FASTQ file(s) into Fastq/")
    for f in fastq_files:
        status = _make_symlink(f.resolve(), fastq_subdir / f.name, overwrite)
        print(f"  [{status}] {f.name}")

    # Symlink the three single files.
    for src, dst in [
        (reference,    ref_subdir  / reference.name),
        (goi,          goi_subdir  / goi.name),
        (sample_sheet, xlsx_subdir / sample_sheet.name),
    ]:
        status = _make_symlink(src.resolve(), dst, overwrite)
        print(f"  [{status}] {dst.parent.name}/{dst.name}")

    print(f"[organize_inputs] Input structure ready at: {run_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build the RNAseq-var input directory structure using symlinks.",
    )
    parser.add_argument("--fastq-dir",    required=True,
                        help="Directory containing paired FASTQ files.")
    parser.add_argument("--sample-sheet", required=True,
                        help="Sample sheet Excel file (.xlsx).")
    parser.add_argument("--reference",    required=True,
                        help="Reference genome FASTA.")
    parser.add_argument("--goi",          required=True,
                        help="Gene-of-interest FASTA.")
    parser.add_argument("--run-dir",      required=True,
                        help="Destination run directory (created if absent).")
    parser.add_argument("--overwrite", action="store_true",
                        help="Replace existing symlinks.")
    args = parser.parse_args()

    run_dir = pathlib.Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    try:
        organize(
            fastq_dir=pathlib.Path(args.fastq_dir).resolve(),
            sample_sheet=pathlib.Path(args.sample_sheet).resolve(),
            reference=pathlib.Path(args.reference).resolve(),
            goi=pathlib.Path(args.goi).resolve(),
            run_dir=run_dir,
            overwrite=args.overwrite,
        )
    except ValueError as exc:
        print(f"[organize_inputs] ERROR: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
