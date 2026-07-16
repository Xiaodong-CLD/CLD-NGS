#!/usr/bin/env python3
"""Allele-account the GATK-uninformative reads at each called site.

Each read in the uninformative BAM carries a ``UV`` tag (``chrom:pos:ref>alt``)
naming the GATK site it was dropped at. For every read we ask one question --
what allele does it show at that site -- and bucket it into exactly three
outcomes using deterministic CIGAR geometry only:

  - ALT          : the read reproduces the called event (deletion gap over the
                   called span, inserted bases at the anchor, or the ALT base).
  - REF          : the site is spanned by plain aligned reference, no event.
  - uncountable  : the site falls in a soft clip / past the read end, or the
                   read carries a different event. This is the artifact sink --
                   it is never subdivided, so sequencing noise cannot force a
                   judgment call.

Outputs one summary row per GATK ALT (UV site). Because 'ALT' means a different
allele at each site, sites are never pooled -- a sample with two GATK calls
yields two rows.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys
from typing import List, Optional, Tuple

import pysam

from pipeline_io import log

EXCLUDE_FLAGS = 4 | 256 | 1024 | 2048  # 3332: unmapped/secondary/dup/supp
CIGAR_S = 4


def _get_softclip_lengths(cigar: List[Tuple[int, int]]) -> Tuple[int, int]:
    """Return (5p_len, 3p_len) from CIGAR tuples."""
    len_5p, len_3p = 0, 0
    if cigar:
        if cigar[0][0] == CIGAR_S:
            len_5p = cigar[0][1]
        if cigar[-1][0] == CIGAR_S:
            len_3p = cigar[-1][1]
    return len_5p, len_3p


def _parse_uv(uv: Optional[str]) -> Optional[Tuple[str, int, str, str]]:
    """Parse a 'chrom:pos:ref>alt' UV tag into (chrom, pos1, ref, alt) or None."""
    if not uv or ">" not in uv:
        return None
    try:
        head, alt = uv.rsplit(">", 1)
        chrom, pos_str, ref = head.split(":", 2)
        return chrom, int(pos_str), ref, alt
    except ValueError:
        return None


def _aligned_base_at(read: pysam.AlignedSegment, pos0: int) -> Optional[str]:
    """Read base aligned to 0-based reference pos0 (M columns only), or None."""
    seq = read.query_sequence
    if seq is None:
        return None
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if rpos == pos0:
            return seq[qpos]
    return None


def _ref_region_aligned(read: pysam.AlignedSegment, start0: int, end0: int) -> bool:
    """True if [start0, end0) is spanned by a single aligned (M) block."""
    for b_start, b_end in read.get_blocks():
        if b_start <= start0 and b_end >= end0:
            return True
    return False


def _has_deletion_overlapping(read: pysam.AlignedSegment, ds: int, de: int) -> bool:
    """True if a CIGAR deletion (D) reference span overlaps [ds, de).

    Walks the CIGAR for D operations directly rather than inferring gaps from
    get_blocks(), because an insertion (I) splits the aligned blocks with a
    zero-width reference gap that would otherwise be mistaken for a deletion.
    """
    rpos = read.reference_start
    if rpos is None:
        return False
    for op, length in (read.cigartuples or []):
        if op == 2:  # D -- deletion spans [rpos, rpos + length)
            if rpos < de and ds < rpos + length:
                return True
            rpos += length
        elif op in (0, 7, 8, 3):  # M/=/X/N consume reference
            rpos += length
        # I/S/H/P consume no reference
    return False


def _has_insertion_at(read: pysam.AlignedSegment, ref_pos0: int, tol: int = 2) -> bool:
    """True if the read has an insertion within tol bp of reference pos ref_pos0."""
    rpos = read.reference_start
    if rpos is None:
        return False
    for op, length in (read.cigartuples or []):
        if op in (0, 7, 8, 2, 3):  # M/=/X/D/N consume reference
            rpos += length
        elif op == 1:  # I -- sits between rpos-1 and rpos
            if abs(rpos - ref_pos0) <= tol:
                return True
    return False


def _ref_subclass(read: pysam.AlignedSegment, clip_bp: int) -> str:
    """Sub-label a REF read: 'REF_clipped' if it carries a long soft clip
    (a cut-site junction read that is reference at the called site but edited
    elsewhere), else 'REF_clean' (genuine wild-type span).
    """
    len_5p, len_3p = _get_softclip_lengths(read.cigartuples or [])
    return "REF_clipped" if max(len_5p, len_3p) >= clip_bp else "REF_clean"


def _uncountable_reason(read: pysam.AlignedSegment, ds: int, de: int) -> str:
    """Explain why a read could not be called at the site [ds, de).

    'not_covered'  : the site lies outside the read's aligned reference span
                     and no soft clip sits on that side (short read / mate).
    'in_softclip'  : the site is past the aligned span on a side that is soft
                     clipped -- the read's bases over the site were clipped.
    'other_indel'  : the site is within the aligned span but broken by a
                     different indel (not the called event).
    """
    rs, re = read.reference_start, read.reference_end
    if rs is None or re is None:
        return "not_covered"
    if ds < rs or de > re:
        ct = read.cigartuples or []
        has_5p_clip = bool(ct) and ct[0][0] == CIGAR_S
        has_3p_clip = bool(ct) and ct[-1][0] == CIGAR_S
        if (de <= rs and has_5p_clip) or (ds >= re and has_3p_clip):
            return "in_softclip"
        return "not_covered"
    return "other_indel"


def classify_allele(
    read: pysam.AlignedSegment, pos1: int, ref: str, alt: str, ref_clip_bp: int = 20,
) -> Tuple[str, str]:
    """Classify a read at its called site.

    Returns (allele, detail). ``allele`` is the unchanged primary verdict
    'ALT' / 'REF' / 'uncountable' (deterministic CIGAR geometry only). ``detail``
    is a descriptive sub-label: for REF it is 'REF_clipped' / 'REF_clean'; for
    uncountable it is the reason ('in_softclip' / 'not_covered' / 'other_indel'
    / 'other_allele'); for ALT it is ''.
    """
    pos0 = pos1 - 1
    if len(ref) > len(alt):            # deletion (left-anchored VCF form)
        ds = pos0 + len(alt)
        de = pos0 + len(ref)
        if _has_deletion_overlapping(read, ds, de):
            return "ALT", ""
        if _ref_region_aligned(read, ds, de):
            return "REF", _ref_subclass(read, ref_clip_bp)
        return "uncountable", _uncountable_reason(read, ds, de)
    if len(alt) > len(ref):            # insertion (left-anchored VCF form)
        anchor_end0 = pos0 + len(ref)
        if _has_insertion_at(read, anchor_end0):
            return "ALT", ""
        if _ref_region_aligned(read, pos0, pos0 + len(ref)):
            return "REF", _ref_subclass(read, ref_clip_bp)
        return "uncountable", _uncountable_reason(read, pos0, pos0 + len(ref))
    # substitution / MNV (equal lengths)
    bases = []
    for k in range(len(ref)):
        b = _aligned_base_at(read, pos0 + k)
        if b is None:
            return "uncountable", _uncountable_reason(read, pos0, pos0 + len(ref))
        bases.append(b.upper())
    obs = "".join(bases)
    if obs == alt.upper():
        return "ALT", ""
    if obs == ref.upper():
        return "REF", _ref_subclass(read, ref_clip_bp)
    return "uncountable", "other_allele"


def main(argv: List[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Allele-account GATK-uninformative reads (ALT/REF/uncountable)."
    )
    ap.add_argument(
        "-b", "--bam", required=True, type=pathlib.Path,
        help="Uninformative BAM from isolate_uninformative_reads.py (UV-tagged)."
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
        help="Max softclip >= this sub-labels a REF read 'REF_clipped' (default: 20)."
    )
    args = ap.parse_args(argv)

    if not args.bam.exists():
        log(f"ERROR: BAM not found: {args.bam}")
        return 2

    sample = args.sample or args.bam.stem.replace(".sorted", "")
    log(f"sample = {sample}")

    # Per-site accumulators keyed by the UV tag string. Each GATK ALT (UV site)
    # gets its own summary row, because 'ALT' means a different allele at each
    # site and must not be pooled across sites.
    def _new_site(chrom, pos1, ref, alt):
        return {
            "chrom": chrom, "pos": pos1, "ref": ref, "alt": alt,
            "ALT": 0, "REF": 0, "uncountable": 0,
            "REF_clipped": 0, "REF_clean": 0,
            "in_softclip": 0, "not_covered": 0, "other_indel": 0, "other_allele": 0,
        }

    site_acc: dict = {}
    n_missing_uv = 0
    n_sample_reads = 0  # all primary reads passing EXCLUDE_FLAGS (sample denominator)

    in_bam = pysam.AlignmentFile(str(args.bam), "rb")
    try:
        for read in in_bam.fetch(until_eof=True):
            if read.flag & EXCLUDE_FLAGS:
                continue
            n_sample_reads += 1
            try:
                uv = read.get_tag("UV")
            except (KeyError, TypeError):
                uv = None
            parsed = _parse_uv(uv) if uv else None
            if parsed is None:
                n_missing_uv += 1
                continue
            _chrom, pos1, ref, alt = parsed
            acc = site_acc.get(uv)
            if acc is None:
                acc = site_acc[uv] = _new_site(_chrom, pos1, ref, alt)
            allele, detail = classify_allele(read, pos1, ref, alt, args.long_softclip_bp)
            acc[allele] += 1
            if allele == "REF":
                acc[detail] += 1
            elif allele == "uncountable":
                acc[detail] += 1
    finally:
        in_bam.close()

    # Roll site accumulators into one summary row per UV site (sorted by pos).
    sample_total = n_sample_reads if n_sample_reads > 0 else 1
    site_rows: List[dict] = []
    tot_alt = tot_ref = 0
    for uv, acc in sorted(site_acc.items(), key=lambda kv: (kv[1]["chrom"], kv[1]["pos"])):
        n_alt, n_ref, n_unc = acc["ALT"], acc["REF"], acc["uncountable"]
        site_total = n_alt + n_ref + n_unc
        callable_denom = n_alt + n_ref
        alt_fraction = n_alt / callable_denom if callable_denom > 0 else 0.0
        tot_alt += n_alt
        tot_ref += n_ref
        site_rows.append({
            "sample": sample,
            "chrom": acc["chrom"],
            "pos": acc["pos"],
            "ref": acc["ref"],
            "alt": acc["alt"],
            "total_uninformative_reads": site_total,
            "sample_read_frac": site_total / sample_total,
            "n_alt": n_alt,
            "n_ref": n_ref,
            "n_uncountable": n_unc,
            "n_ref_clipped": acc["REF_clipped"],
            "n_ref_clean": acc["REF_clean"],
            "uninf_alt_frac": alt_fraction,
        })

    grand_total = tot_alt + tot_ref + sum(a["uncountable"] for a in site_acc.values())
    log(f"UV sites = {len(site_rows)}; total reads (post-filter) = {grand_total}"
        f"{f'; reads missing UV tag = {n_missing_uv}' if n_missing_uv else ''}")
    for row in site_rows:
        log(f"  {row['chrom']}:{row['pos']} {row['ref']}>{row['alt']}: "
            f"ALT={row['n_alt']} REF={row['n_ref']} "
            f"(clipped={row['n_ref_clipped']}/clean={row['n_ref_clean']}) "
            f"uncountable={row['n_uncountable']} alt_frac={row['uninf_alt_frac']:.4f}")

    fieldnames = [
        "sample", "chrom", "pos", "ref", "alt",
        "total_uninformative_reads", "sample_read_frac",
        "n_alt", "n_ref", "n_uncountable",
        "n_ref_clipped", "n_ref_clean",
        "uninf_alt_frac",
    ]

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in site_rows:
            w.writerow(row)
    log(f"Wrote {args.output_tsv} ({len(site_rows)} site rows)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
