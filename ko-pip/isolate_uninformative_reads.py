#!/usr/bin/env python3
"""Isolate GATK-uninformative reads for IGV inspection.

For every PASS / '.' variant in a GATK VCF, emit the reads that overlap
the variant position in the raw sorted BAM but are absent from the GATK
HaplotypeCaller -bamout at that site (i.e. covered the locus but did not
contribute to any AD bucket). Output:

  <output_dir>/uninformative.sorted.bam   (+ .bai)
  <output_dir>/uninformative_sites.tsv    (one row per (read,site))

The companion notebook (expt/uninformative_read_analysis/evaluate.ipynb)
joins the TSV with per-read CIGAR shape to bucket the reasons.

Read-universe filter: SAM flag -F 3332 ->
exclude unmapped (4), secondary (256), duplicate (1024), supplementary
(2048). Read group / MAPQ are kept as-is so we can see HC-MAPQ exclusions.

Coordinates: external (VCF, TSV) are 1-based inclusive; internal pysam
indices are 0-based half-open. Conversion only at I/O.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import subprocess
import sys
import tempfile
from typing import Dict, List, Set, Tuple

import pysam

LOG_PREFIX = "[Takeda Pipeline Message]"

# SAM flags to exclude (`samtools view -F 3332`).
EXCLUDE_FLAGS = 4 | 256 | 1024 | 2048  # = 3332


def log(msg: str) -> None:
    print(f"{LOG_PREFIX} {msg}", flush=True)


def _read_passes_filter(read: pysam.AlignedSegment) -> bool:
    """Mirror `samtools view -F 3332`."""
    return (read.flag & EXCLUDE_FLAGS) == 0


def _ensure_bamout(
    bam_path: pathlib.Path,
    reference: pathlib.Path,
    output_dir: pathlib.Path,
    bamout_arg: pathlib.Path | None,
    sample: str,
    threads: int,
) -> pathlib.Path:
    """Return path to an indexed HC bamout; regenerate if missing.

    VCF output is written to a temp file and discarded after indexing.
    """
    if bamout_arg is not None and bamout_arg.exists():
        bamout = bamout_arg
    else:
        bamout = output_dir / f"{sample}_bamout.bam"
    if not bamout.exists():
        log(f"generating HC bamout -> {bamout}")
        bamout.parent.mkdir(parents=True, exist_ok=True)
        # Write VCF to temp location and discard it after (not needed for analysis)
        with tempfile.TemporaryDirectory() as tmpdir:
            vcf_out = pathlib.Path(tmpdir) / "temp.vcf.gz"
            cmd = [
                "gatk", "HaplotypeCaller",
                "-R", str(reference),
                "-I", str(bam_path),
                "-O", str(vcf_out),
                "-bamout", str(bamout),
                "--max-reads-per-alignment-start", "0",
                "--stand-call-conf", "20",
                "--active-probability-threshold", "0.0002",
                "--min-base-quality-score", "5",
            ]
            subprocess.run(cmd, check=True)
            # VCF is automatically cleaned up when tmpdir exits
    if not (bamout.with_suffix(bamout.suffix + ".bai").exists()
            or bamout.with_suffix(".bai").exists()):
        log(f"indexing bamout {bamout}")
        subprocess.run(["samtools", "index", "-@", str(threads), str(bamout)], check=True)
    return bamout


def _iter_pass_sites(vcf_path: pathlib.Path):
    """Yield (chrom, pos_1based, ref, alt_str) for each PASS / '.' ALT."""
    with pysam.VariantFile(str(vcf_path)) as vf:
        for rec in vf:
            filt_keys = list(rec.filter.keys())
            # Empty FILTER (VCF '.') -> pysam returns []; treat as PASS-equivalent.
            is_pass_or_dot = (not filt_keys) or ("PASS" in filt_keys)
            if not is_pass_or_dot:
                continue
            for alt in (rec.alts or ()):
                yield rec.chrom, rec.pos, rec.ref, alt


def _site_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{chrom}:{pos}:{ref}>{alt}"


def _informative_qnames_at(bamout: pysam.AlignmentFile, chrom: str, pos_1: int) -> Set[Tuple[str, bool]]:
    """qname/is_read1 of reads present in bamout overlapping pos_1.

    HC bamout preserves original qnames for reads it kept; synthetic
    haplotype reads (qname like 'HC*') are also returned but never
    intersect the raw-BAM qname set.
    """
    out: Set[Tuple[str, bool]] = set()
    for r in bamout.fetch(chrom, pos_1 - 1, pos_1):
        if r.qname is None:
            continue
        out.add((r.qname, bool(r.is_read1)))
    return out


def _aligned_length(read: pysam.AlignedSegment) -> int:
    if read.reference_end is None or read.reference_start is None:
        return 0
    return read.reference_end - read.reference_start


def main(argv: List[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Isolate GATK-uninformative reads.")
    p.add_argument("-b", "--bam", required=True, type=pathlib.Path,
                   help="Raw sorted, indexed BAM (e.g. tmp/<sample>/<sample>.sorted.bam)")
    p.add_argument("-v", "--vcf", required=True, type=pathlib.Path,
                   help="GATK HaplotypeCaller VCF (e.g. Results/<sample>/gatk/gatk_variants.vcf.gz)")
    p.add_argument("-r", "--reference", required=True, type=pathlib.Path,
                   help="Reference FASTA (needed only if bamout must be regenerated)")
    p.add_argument("--bamout", type=pathlib.Path, default=None,
                   help="Existing HC -bamout BAM; regenerated next to output if absent.")
    p.add_argument("-o", "--output-dir", required=True, type=pathlib.Path,
                   help="Output directory (created if missing).")
    p.add_argument("--sample", default=None,
                   help="Sample id; defaults to the input BAM stem with .sorted stripped.")
    p.add_argument("--threads", type=int, default=4)
    args = p.parse_args(argv)

    for path in (args.bam, args.vcf):
        if not path.exists():
            print(f"ERROR: missing input {path}", file=sys.stderr)
            return 2
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sample = args.sample or args.bam.stem.replace(".sorted", "")
    log(f"sample = {sample}")

    bamout_path = _ensure_bamout(
        args.bam, args.reference, args.output_dir, args.bamout, sample, args.threads
    )

    # --- Pass 1: per site, find uninformative (qname, is_read1) ----------
    sites = list(_iter_pass_sites(args.vcf))
    log(f"PASS/. variants (per-ALT) = {len(sites)}")
    if not sites:
        log("no PASS variants; nothing to do")
        return 0

    # (qname, is_read1) -> list of site_ids it was uninformative at.
    uninf_to_sites: Dict[Tuple[str, bool], List[str]] = {}
    # Stderr cardinality table: site_id -> (n_raw, n_bamout, n_uninf)
    site_card: List[Tuple[str, int, int, int]] = []

    raw_bam = pysam.AlignmentFile(str(args.bam), "rb")
    bamout = pysam.AlignmentFile(str(bamout_path), "rb")
    try:
        for chrom, pos_1, ref, alt in sites:
            sid = _site_id(chrom, pos_1, ref, alt)
            informative = _informative_qnames_at(bamout, chrom, pos_1)
            n_raw = 0
            n_uninf = 0
            for r in raw_bam.fetch(chrom, pos_1 - 1, pos_1):
                if not _read_passes_filter(r):
                    continue
                n_raw += 1
                key = (r.qname, bool(r.is_read1))
                if key in informative:
                    continue
                n_uninf += 1
                uninf_to_sites.setdefault(key, []).append(sid)
            site_card.append((sid, n_raw, len(informative), n_uninf))
    finally:
        bamout.close()
        raw_bam.close()

    log(f"unique uninformative (qname,mate) pairs = {len(uninf_to_sites)}")
    for sid, nr, nb, nu in site_card:
        print(f"  {sid}\traw={nr}\tbamout={nb}\tuninf={nu}", file=sys.stderr)

    # --- Pass 2: emit BAM + TSV ------------------------------------------
    unsorted_bam = args.output_dir / f"{sample}.bam"
    sorted_bam = args.output_dir / f"{sample}.sorted.bam"

    written_keys: Set[Tuple[str, bool]] = set()

    # Reopen raw BAM for a fresh linear scan; mixing random-access fetch with
    # until_eof on the same handle leaves the iterator at EOF.
    raw_bam = pysam.AlignmentFile(str(args.bam), "rb")
    with pysam.AlignmentFile(str(unsorted_bam), "wb", template=raw_bam) as out_bam:
        for r in raw_bam.fetch(until_eof=True):
            if not _read_passes_filter(r):
                continue
            key = (r.qname, bool(r.is_read1))
            if key not in uninf_to_sites:
                continue
            sids = uninf_to_sites[key]
            first_sid = sids[0]
            cigar_str = r.cigarstring or ""
            aln_len = _aligned_length(r)
            ref_start_1 = (r.reference_start + 1) if r.reference_start is not None else 0
            if key not in written_keys:
                # Tag and write the BAM record exactly once.
                try:
                    r.set_tag("UV", first_sid, value_type="Z")
                    r.set_tag("UR", "not_in_bamout", value_type="Z")
                except Exception:  # extremely defensive; bad tags shouldn't kill the run
                    pass
                out_bam.write(r)
                written_keys.add(key)

    raw_bam.close()

    log(f"sorting + indexing -> {sorted_bam}")
    subprocess.run(
        ["samtools", "sort", "-@", str(args.threads), "-o", str(sorted_bam), str(unsorted_bam)],
        check=True,
    )
    subprocess.run(["samtools", "index", "-@", str(args.threads), str(sorted_bam)], check=True)
    unsorted_bam.unlink(missing_ok=True)

    log(f"BAM records written = {len(written_keys)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
