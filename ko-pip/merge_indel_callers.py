#!/usr/bin/env python3
"""Merge GATK + VarScan2 indel calls into unified table.

Normalize both VCFs to canonical left-aligned form, deduplicate, and combine
evidence into a single TSV.

Coordinates:
  - Output is 1-based inclusive (external format).
  - Internal use is 0-based half-open (pysam convention).

Normalization:
  - Uses bcftools norm -f ref.fa to left-align both VCFs.
  - Extracts INDEL records only (filters out SNPs).

Merging:
  - Canonical key: (chrom, pos, ref, alt)
  - For each unique key, populates GATK columns (if found), VarScan2 columns (if found).
  - Labels source as: gatk_only, varscan2_only, both.

CLI
---
    python merge_indel_callers.py \\
        --gatk-vcf Results/<s>/indel/gatk_variants.vcf.gz \\
        --varscan-indel-vcf Results/<s>/indel/*.indel.vcf \\
        --reference-path data/reference/ref.fa \\
        --chrom Fut8 \\
        --sample-id <s> \\
        --output-dir Results/<s>/indel
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import os
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional, Tuple

import pysam

logger = logging.getLogger("merge_indel_callers")
MIN_TOTAL_INFORMATIVE_FRAC = 0.01  # gate: VarScan depth must be >= this fraction of total mapped reads

MERGED_COLUMNS = [
    "chrom",
    "start_pos",
    "end_pos",
    "seq",
    "type",
    "length",
    "source",
    "gatk_alt_frequency_direct",
    "gatk_sum_ad_frac",
    "varscan2_alt_frequency",
    "varscan2_read_support",
    "varscan2_sdp",
    "varscan2_sum_ad_frac",
    "edge_artifact",
]

# An indel whose 1-based start position falls within this many bp of either
# amplicon terminus is flagged edge_artifact='yes'. These terminal calls (e.g.
# an INS right at the last base) are alignment/primer-edge noise; the flag is
# kept in the per-sample merged_indel.tsv but the run-level aggregator drops
# them from all_indel.tsv.
EDGE_ARTIFACT_MARGIN_BP = 3


def _normalize_vcf(vcf_path: str, reference_path: str, output_path: str) -> None:
    """Left-align and split multiallelics with bcftools norm.

    ``-m-any`` splits multiallelic records into biallelic ones so a GATK
    multiallelic indel (e.g. ``CTTTG>*,C``) and the matching VarScan call
    (``CTTTG>C``) dedup on the same (chrom, pos, ref, alt) key.

    VarScan ``--output-vcf`` files omit ``##contig`` header lines, which makes
    bcftools refuse to parse them. When the reference ``.fai`` is available and
    the VCF has no contig lines, reheader it first to inject them.
    """
    if not os.path.exists(vcf_path):
        logger.warning(f"VCF not found: {vcf_path}")
        return
    try:
        norm_input = vcf_path
        fai_path = reference_path + ".fai"
        if os.path.exists(fai_path) and not _has_contig_header(vcf_path):
            reheadered = output_path + ".reheader.vcf"
            subprocess.run(
                ["bcftools", "reheader", "--fai", fai_path, "-o", reheadered, vcf_path],
                check=True, capture_output=True,
            )
            norm_input = reheadered
        cmd = [
            "bcftools", "norm",
            "-f", reference_path,
            "-m-any",
            "-o", output_path,
            "-O", "z",
            norm_input,
        ]
        subprocess.run(cmd, check=True, capture_output=True)
        # Index the normalized VCF
        subprocess.run(["bcftools", "index", output_path], check=True, capture_output=True)
        logger.info(f"Normalized: {vcf_path} -> {output_path}")
    except subprocess.CalledProcessError as e:
        logger.error(f"bcftools norm failed: {e.stderr.decode()}")
        raise


def _has_contig_header(vcf_path: str) -> bool:
    """True if the VCF declares at least one ``##contig`` header line."""
    try:
        opener = gzip.open if vcf_path.endswith(".gz") else open
        with opener(vcf_path, "rt") as fh:
            for line in fh:
                if not line.startswith("##"):
                    return False  # past the header
                if line.startswith("##contig"):
                    return True
    except OSError:
        return False
    return False


def _is_indel(ref: str, alt: str) -> bool:
    """Check if variant is an indel (not a SNP)."""
    return len(ref) != len(alt)


def _pct_to_fraction(value: str) -> float:
    """Parse a VarScan VarFreq string like '95.84%' into a 0..1 fraction."""
    s = (value or "").strip()
    if s.endswith("%"):
        s = s[:-1]
    try:
        return float(s) / 100.0
    except ValueError:
        return 0.0


def _fmt_float(value, ndigits: int) -> str:
    """Format a numeric value to fixed precision; empty string if not numeric."""
    if isinstance(value, (int, float)):
        return f"{value:.{ndigits}f}"
    return "."


def _parse_gatk_vcf(vcf_path: str) -> Dict[Tuple[str, int, str, str], Dict]:
    """Extract indels from GATK VCF. Return dict keyed by (chrom, pos, ref, alt)."""
    gatk_indels = {}
    if not os.path.exists(vcf_path):
        return gatk_indels
    try:
        vcf = pysam.VariantFile(vcf_path)
        for rec in vcf:
            # Only process indels on target contig
            ref = rec.ref
            alt = rec.alts[0] if rec.alts else ref
            if alt == "*":
                continue  # spanning-deletion placeholder from multiallelic split
            if not _is_indel(ref, alt):
                continue
            
            # Extract fields
            ac = rec.info.get("AC", [0])
            ac = ac[0] if isinstance(ac, (list, tuple)) else ac
            af = rec.info.get("AF", [0.0])
            af = af[0] if isinstance(af, (list, tuple)) else af
            dp = rec.info.get("DP", 0) if "DP" in rec.info else 0
            # AD is a per-sample FORMAT field; read it from the sample, not the
            # header metadata (rec.format is header metadata and not subscriptable).
            ad = (0, 0)
            if len(rec.samples) > 0:
                sample = rec.samples[list(rec.samples)[0]]
                ad_val = sample.get("AD")
                if ad_val:
                    ad = ad_val
            
            # Compute alt_depth and alt_frequency from AD
            alt_depth = ad[1] if len(ad) > 1 else 0
            total_ad = sum(ad) if ad else 0
            alt_freq_direct = (alt_depth / total_ad) if total_ad else 0.0
            # pct_sum_ad: fraction of raw locus depth (INFO/DP) GATK attributed to
            # an allele (sum_AD / DP_INFO). This is the GATK keep-gate metric.
            pct_sum_ad = (total_ad / dp) if dp else 0.0

            qual = rec.qual if rec.qual else 0.0

            key = (rec.contig, rec.start, ref, alt)  # rec.start is 0-based
            gatk_indels[key] = {
                "allele_count": ac,
                "allele_frequency": af,
                "alt_depth_sum": alt_depth,
                "alt_frequency_direct": alt_freq_direct,
                "pileup_dp": dp,
                "qual": qual,
                "sum_ad": total_ad,
                "pct_sum_ad": pct_sum_ad,
            }
        vcf.close()
    except Exception as e:
        logger.error(f"Failed to parse GATK VCF {vcf_path}: {e}")
    return gatk_indels


def _parse_varscan_vcf(vcf_path: str) -> Dict[Tuple[str, int, str, str], Dict]:
    """Extract indels from a VarScan2 ``--output-vcf`` indel VCF.

    VarScan stores per-sample depths and frequency in the FORMAT/sample
    column (not INFO): ``RD`` = ref depth, ``AD`` = alt depth (a single
    scalar, unlike GATK's [ref, alt] pair), ``FREQ`` = a percent string such
    as ``"43.15%"``. Keys are (chrom, rec.start, ref, alt) with rec.start
    0-based, matching the GATK parser so the two callers dedup.
    """
    varscan_indels: Dict[Tuple[str, int, str, str], Dict] = {}
    if not vcf_path or not os.path.exists(vcf_path):
        return varscan_indels
    try:
        vcf = pysam.VariantFile(vcf_path)
        for rec in vcf:
            ref = rec.ref
            alt = rec.alts[0] if rec.alts else ref
            if alt == "*":
                continue
            if not _is_indel(ref, alt):
                continue

            freq_val = 0.0
            rd = 0
            ad = 0
            sdp = 0
            if len(rec.samples) > 0:
                sample = rec.samples[list(rec.samples)[0]]
                freq = sample.get("FREQ")
                if isinstance(freq, (list, tuple)):
                    freq = freq[0] if freq else ""
                freq_val = _pct_to_fraction(freq if isinstance(freq, str) else "")
                rd_val = sample.get("RD")
                rd = rd_val[0] if isinstance(rd_val, (list, tuple)) and rd_val else (rd_val or 0)
                ad_val = sample.get("AD")
                ad = ad_val[0] if isinstance(ad_val, (list, tuple)) and ad_val else (ad_val or 0)
                sdp_val = sample.get("SDP")
                sdp = sdp_val[0] if isinstance(sdp_val, (list, tuple)) and sdp_val else (sdp_val or 0)

            # pct_sum_ad: fraction of raw read depth (SDP) attributed to ref+alt
            # ((AD + RD) / SDP). VarScan mirror of the GATK keep-gate metric.
            attributed = (ad or 0) + (rd or 0)
            pct_sum_ad = (attributed / sdp) if sdp else 0.0

            key = (rec.contig, rec.start, ref, alt)  # rec.start is 0-based
            varscan_indels[key] = {
                "alt_frequency": freq_val,
                "read_support": f"{ad}:{rd}",  # alt:ref
                "sdp": sdp,
                "pct_sum_ad": pct_sum_ad,
                "attributed": attributed,  # ad + rd; used for total-reads gate
            }
        vcf.close()
    except Exception as e:
        logger.error(f"Failed to parse VarScan2 VCF {vcf_path}: {e}")
    return varscan_indels


def _indel_type(ref: str, alt: str) -> str:
    """Return DEL or INS."""
    return "DEL" if len(ref) > len(alt) else "INS"


def _indel_length(ref: str, alt: str) -> int:
    """Return absolute indel length."""
    return abs(len(ref) - len(alt))


def _is_edge_artifact(pos_1based: int, ref_len: int, margin: int = EDGE_ARTIFACT_MARGIN_BP) -> bool:
    """True if a 1-based indel start sits within ``margin`` bp of either end.

    With margin=3 on a 130 bp amplicon this flags positions 1-3 (5' end) and
    128-130 (3' end). When ref_len is unknown (0) nothing is flagged.
    """
    if not ref_len:
        return False
    return pos_1based <= margin or pos_1based > ref_len - margin


def merge_indel_callers(
    gatk_vcf: Optional[str],
    varscan_indel_vcf: Optional[str],
    reference_path: str,
    chrom: str,
    sample_id: str,
    output_dir: str,
    total_mapped_reads: Optional[int] = None,
) -> str:
    """Merge GATK and VarScan2 indels into unified TSV.

    ``total_mapped_reads`` (mapped reads in the sorted BAM) is recorded only as
    a ``# total_mapped_reads`` comment banner for the per-run aggregator's
    per-sample header; it does not affect any row.

    Returns path to merged TSV.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    out_path = os.path.join(output_dir, "merged_indel.tsv")

    # GATK is a real VCF: left-align + split multiallelics, then parse.
    # VarScan2 (--output-vcf) is also a real VCF, so normalize it the same
    # way; both then dedup on a shared left-aligned (chrom, pos, ref, alt) key.
    with tempfile.TemporaryDirectory() as tmpdir:
        if gatk_vcf and os.path.exists(gatk_vcf):
            gatk_norm = os.path.join(tmpdir, "gatk_norm.vcf.gz")
            _normalize_vcf(gatk_vcf, reference_path, gatk_norm)
            gatk_indels = _parse_gatk_vcf(gatk_norm)
        else:
            gatk_indels = {}

        if varscan_indel_vcf and os.path.exists(varscan_indel_vcf):
            varscan_norm = os.path.join(tmpdir, "varscan_norm.vcf.gz")
            _normalize_vcf(varscan_indel_vcf, reference_path, varscan_norm)
            varscan_indels = _parse_varscan_vcf(varscan_norm)
        else:
            varscan_indels = {}

    # Contig lengths for edge-artifact flagging (position near amplicon end).
    ref_lengths: Dict[str, int] = {}
    try:
        with pysam.FastaFile(reference_path) as _fa:
            ref_lengths = dict(zip(_fa.references, _fa.lengths))
    except Exception:
        ref_lengths = {}

    # Merge
    all_keys = set(gatk_indels.keys()) | set(varscan_indels.keys())
    
    out_rows: List[Dict] = []
    for key in sorted(all_keys):
        chrom_key, pos_0based, ref, alt = key
        
        # Convert to 1-based for output
        pos_1based = pos_0based + 1
        
        gatk_data = gatk_indels.get(key, {})
        varscan_data = varscan_indels.get(key, {})

        # Drop VarScan evidence where informative depth < 1% of total mapped reads.
        if varscan_data and total_mapped_reads and \
                varscan_data.get("attributed", 0) < MIN_TOTAL_INFORMATIVE_FRAC * total_mapped_reads:
            varscan_data = {}

        in_gatk = bool(gatk_data)
        in_varscan = bool(varscan_data)

        if not in_gatk and not in_varscan:
            continue  # VarScan-only row dropped by the depth gate above

        if in_gatk and in_varscan:
            source = "both"
        elif in_gatk:
            source = "gatk_only"
        else:
            source = "varscan2_only"
        
        vtype = _indel_type(ref, alt)
        if vtype == "DEL":
            start_pos = pos_1based + 1
            end_pos = pos_1based + len(ref) - 1
            seq = ref[1:]
        else:  # INS
            start_pos = pos_1based
            end_pos = pos_1based
            seq = alt[1:]

        row = {
            "chrom": chrom_key,
            "start_pos": start_pos,
            "end_pos": end_pos,
            "seq": seq,
            "type": vtype,
            "length": _indel_length(ref, alt),
            "source": source,
            "gatk_alt_frequency_direct": _fmt_float(gatk_data.get("alt_frequency_direct"), 6),
            "gatk_sum_ad_frac": _fmt_float(gatk_data.get("pct_sum_ad"), 6),
            "varscan2_alt_frequency": _fmt_float(varscan_data.get("alt_frequency"), 6),
            "varscan2_read_support": varscan_data.get("read_support", "."),
            "varscan2_sdp": varscan_data.get("sdp", "."),
            "varscan2_sum_ad_frac": _fmt_float(varscan_data.get("pct_sum_ad"), 6),
            "edge_artifact": "yes" if _is_edge_artifact(
                pos_1based, ref_lengths.get(chrom_key, 0)
            ) else "no",
        }
        out_rows.append(row)
    
    # Write output
    with open(out_path, "w", newline="") as fh:
        fh.write(f"# sample: {sample_id}\n")
        if total_mapped_reads is not None:
            fh.write(f"# total_mapped_reads: {total_mapped_reads}\n")
        w = csv.DictWriter(fh, fieldnames=MERGED_COLUMNS, delimiter="\t")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    
    logger.info(f"Merged {len(gatk_indels)} GATK indels + {len(varscan_indels)} VarScan2 indels -> {len(out_rows)} unique events: {out_path}")
    return out_path


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gatk-vcf", required=False, help="GATK VCF (optional)")
    p.add_argument("--varscan-indel-vcf", required=False, help="VarScan2 indel VCF (optional)")
    p.add_argument("--reference-path", required=True, help="Reference FASTA")
    p.add_argument("--chrom", required=True, help="Target contig")
    p.add_argument("--sample-id", required=True)
    p.add_argument("--output-dir", required=True)
    return p


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)
    out = merge_indel_callers(
        gatk_vcf=args.gatk_vcf,
        varscan_indel_vcf=args.varscan_indel_vcf,
        reference_path=args.reference_path,
        chrom=args.chrom,
        sample_id=args.sample_id,
        output_dir=args.output_dir,
    )
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
