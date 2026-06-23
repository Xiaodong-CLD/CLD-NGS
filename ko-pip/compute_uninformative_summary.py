#!/usr/bin/env python3
"""Summarize uninformative reads by softclip length and mismatch rate.

Reads the uninformative BAM and categorizes each unique read by:
  - max_softclip_len: maximum of 5' and 3' softclip lengths
  - mismatch_rate: NM tag / aligned_length (if available)

Categories:
  - long_softclip: max softclip >= threshold (default 20 bp)
  - high_mismatch_low_softclip: high mismatch rate AND short softclips
  - other: remaining

Outputs one summary row per sample with counts and percentages to TSV.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys
from typing import List, Optional, Tuple

import pysam

LOG_PREFIX = "[Takeda Pipeline Message]"
EXCLUDE_FLAGS = 4 | 256 | 1024 | 2048  # 3332: unmapped/secondary/dup/supp
CIGAR_S = 4


def log(msg: str) -> None:
    print(f"{LOG_PREFIX} {msg}", flush=True)


def _get_softclip_lengths(cigar: List[Tuple[int, int]]) -> Tuple[int, int]:
    """Return (5p_len, 3p_len) from CIGAR tuples."""
    len_5p, len_3p = 0, 0
    if cigar:
        if cigar[0][0] == CIGAR_S:
            len_5p = cigar[0][1]
        if cigar[-1][0] == CIGAR_S:
            len_3p = cigar[-1][1]
    return len_5p, len_3p


def _get_aligned_length(read: pysam.AlignedSegment) -> int:
    """Return alignment span in reference."""
    if read.reference_end is None or read.reference_start is None:
        return 0
    return read.reference_end - read.reference_start


def _categorize_read(
    read: pysam.AlignedSegment,
    long_softclip_bp: int,
    little_softclip_bp: int,
    high_mismatch_rate: float,
) -> Tuple[str, int, float, int]:
    """Categorize one read by softclip and mismatch.

    Returns (category, max_softclip_len, mismatch_rate, aligned_len).
    category in ["long_softclip", "high_mismatch_low_softclip", "other"]
    mismatch_rate = -1.0 if not evaluable.
    """
    if read.cigartuples is None or read.query_length is None:
        return "other", 0, -1.0, 0

    len_5p, len_3p = _get_softclip_lengths(read.cigartuples)
    max_sc = max(len_5p, len_3p)
    aln_len = _get_aligned_length(read)

    # Try to compute mismatch rate from NM tag
    mismatch_rate = -1.0
    try:
        nm = read.get_tag("NM")
        if aln_len > 0:
            mismatch_rate = nm / aln_len
    except (KeyError, TypeError):
        pass

    # Long softclip category
    if max_sc >= long_softclip_bp:
        return "long_softclip", max_sc, mismatch_rate, aln_len

    # High mismatch with little softclip
    if max_sc <= little_softclip_bp and mismatch_rate >= 0 and mismatch_rate >= high_mismatch_rate:
        return "high_mismatch_low_softclip", max_sc, mismatch_rate, aln_len

    return "other", max_sc, mismatch_rate, aln_len


def main(argv: List[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Summarize uninformative reads by softclip and mismatch."
    )
    ap.add_argument(
        "-b", "--bam", required=True, type=pathlib.Path,
        help="Uninformative BAM from isolate_uninformative_reads.py"
    )
    ap.add_argument(
        "-o", "--output-tsv", required=True, type=pathlib.Path,
        help="Output summary TSV (one row)."
    )
    ap.add_argument(
        "--sample", default=None,
        help="Sample id; defaults to BAM stem with .sorted stripped."
    )
    ap.add_argument(
        "--long-softclip-bp", type=int, default=20,
        help="Threshold for long softclip (default: 20)."
    )
    ap.add_argument(
        "--little-softclip-bp", type=int, default=5,
        help="Threshold for little softclip (default: 5)."
    )
    ap.add_argument(
        "--high-mismatch-rate", type=float, default=0.05,
        help="Threshold for high mismatch rate (default: 0.05 = 5%)."
    )
    args = ap.parse_args(argv)

    if not args.bam.exists():
        log(f"ERROR: BAM not found: {args.bam}")
        return 2

    sample = args.sample or args.bam.stem.replace(".sorted", "")
    log(f"sample = {sample}")

    # Count categories
    categories = {
        "long_softclip": 0,
        "high_mismatch_low_softclip": 0,
        "other": 0,
    }
    total_reads = 0
    n_missing_nm = 0

    with pysam.AlignmentFile(str(args.bam), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.flag & EXCLUDE_FLAGS:
                continue
            total_reads += 1
            cat, _, mr, _ = _categorize_read(
                read,
                args.long_softclip_bp,
                args.little_softclip_bp,
                args.high_mismatch_rate,
            )
            categories[cat] += 1
            if mr < 0:
                n_missing_nm += 1

    log(f"Total uninformative reads (post-filter) = {total_reads}")
    log(f"  long_softclip = {categories['long_softclip']}")
    log(f"  high_mismatch_low_softclip = {categories['high_mismatch_low_softclip']}")
    log(f"  other = {categories['other']}")
    log(f"  (reads missing NM tag = {n_missing_nm})")

    # Write summary
    denom = total_reads if total_reads > 0 else 1
    summary_row = {
        "sample": sample,
        "total_uninformative_reads": total_reads,
        "n_long_softclip": categories["long_softclip"],
        "pct_long_softclip": 100.0 * categories["long_softclip"] / denom,
        "n_high_mismatch_low_softclip": categories["high_mismatch_low_softclip"],
        "pct_high_mismatch_low_softclip": 100.0 * categories["high_mismatch_low_softclip"] / denom,
        "n_other": categories["other"],
        "pct_other": 100.0 * categories["other"] / denom,
        "threshold_long_softclip_bp": args.long_softclip_bp,
        "threshold_little_softclip_bp": args.little_softclip_bp,
        "threshold_high_mismatch_rate": args.high_mismatch_rate,
        "n_missing_nm_tag": n_missing_nm,
    }

    fieldnames = [
        "sample",
        "total_uninformative_reads",
        "n_long_softclip",
        "pct_long_softclip",
        "n_high_mismatch_low_softclip",
        "pct_high_mismatch_low_softclip",
        "n_other",
        "pct_other",
        "threshold_long_softclip_bp",
        "threshold_little_softclip_bp",
        "threshold_high_mismatch_rate",
        "n_missing_nm_tag",
    ]

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerow(summary_row)

    log(f"Wrote {args.output_tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
