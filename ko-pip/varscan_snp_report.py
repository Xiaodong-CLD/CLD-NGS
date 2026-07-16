"""Thin VarScan SNP report.

Entry point:
    * ``build_snp_report()`` — emits ``<output_dir>/varscan_snp_report.tsv`` from
        the VarScan ``mpileup2snp`` SNP table. The file is a direct parsed summary,
        not a verdict layer: it keeps the basic support and frequency fields and
        leaves thresholding to downstream aggregation or human review.
"""

from __future__ import annotations

import argparse
import csv
import os
from typing import Iterable

SNP_REPORT_COLUMNS = [
    "sample_id",
    "chrom",
    "position",
    "ref",
    "alt",
    "consensus",
    "ref_support",
    "alt_support",
    "total_informative",
    "var_freq",
    "p_value",
    "alt_strands_plus",
    "alt_strands_minus",
]

OUTPUT_FILE = "varscan_snp_report.tsv"
# A variant site must have at least this fraction of total mapped reads as
# informative (ref + alt) depth to be reported. Sites below this are likely
# low-coverage artefacts (primer ends, large-deletion gaps) rather than real
# variants. Gate is skipped when total_reads is not provided (== 0).
MIN_TOTAL_INFORMATIVE_FRAC = 0.01

def _try_int(x: str) -> int:
    try:
        return int(x)
    except (TypeError, ValueError):
        return 0


def _try_float(x: str) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def _parse_var_freq(token: str) -> float:
    """VarScan emits VarFreq as '12.34%'. Return a fraction in [0,1]."""
    s = (token or "").strip()
    if not s:
        return 0.0
    if s.endswith("%"):
        return _try_float(s[:-1]) / 100.0
    return _try_float(s)


def _parse_cons_reads_freq_pvalue(packed: str) -> tuple[str, int, int, float, str]:
    """Parse VarScan packed field Cons:Cov:Reads1:Reads2:Freq:P-value."""
    items = (packed or "").split(":")
    if len(items) < 6:
        return "", 0, 0, 0.0, ""
    cons = items[0]
    reads1 = _try_int(items[2])
    reads2 = _try_int(items[3])
    var_freq = _parse_var_freq(items[4])
    p_value = items[5]
    return cons, reads1, reads2, var_freq, p_value


def _parse_strand_filter(packed: str) -> tuple[str, str]:
    """Parse VarScan packed field StrandFilter:R1+:R1-:R2+:R2-:pval."""
    items = (packed or "").split(":")
    if len(items) < 5:
        return "", ""
    return items[3], items[4]


def _read_varscan_rows(path: str) -> list[dict]:
    """Read the VarScan SNP tab table. Returns [] if file is missing."""
    if not os.path.exists(path):
        return []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def _build_row(sample_id: str, vs: dict) -> dict:
    chrom = vs.get("Chrom", "")

    # Support both schemas:
    # 1) flat/annotated fields (Reads1/Reads2/VarFreq/...)
    # 2) packed VarScan mpileup2snp fields (Cons:Cov:Reads1:Reads2:Freq:P-value, ...)
    ref_support = _try_int(vs.get("Reads1", "0"))
    alt_support = _try_int(vs.get("Reads2", "0"))
    var_freq = _parse_var_freq(vs.get("VarFreq", ""))
    consensus = vs.get("Cons", "")
    p_value = vs.get("Pvalue", "")
    alt_strands_plus = vs.get("Reads2Plus", "")
    alt_strands_minus = vs.get("Reads2Minus", "")

    if not any([ref_support, alt_support]) and var_freq == 0.0 and not consensus:
        consensus, ref_support, alt_support, var_freq, p_value = _parse_cons_reads_freq_pvalue(
            vs.get("Cons:Cov:Reads1:Reads2:Freq:P-value", "")
        )
    if not alt_strands_plus and not alt_strands_minus:
        alt_strands_plus, alt_strands_minus = _parse_strand_filter(
            vs.get("StrandFilter:R1+:R1-:R2+:R2-:pval", "")
        )

    total = ref_support + alt_support

    return {
        "sample_id": sample_id,
        "chrom": chrom,
        "position": vs.get("Position", ""),
        "ref": vs.get("Ref", ""),
        "alt": vs.get("VarAllele", "") or vs.get("Var", ""),
        "consensus": consensus,
        "ref_support": ref_support,
        "alt_support": alt_support,
        "total_informative": total,
        "var_freq": f"{var_freq:.4f}",
        "p_value": p_value,
        "alt_strands_plus": alt_strands_plus,
        "alt_strands_minus": alt_strands_minus,
    }


def _pick_input(varscan_dir: str, sample_id: str) -> tuple[str, str]:
    """Prefer the amino-annotated file; fall back to the plain SNP table.

    Returns (path, source_label).
    """
    annotated = os.path.join(varscan_dir, f"{sample_id}.snp.amino_annotated.vcf")
    plain = os.path.join(varscan_dir, f"{sample_id}.snp.vcf")
    if os.path.exists(annotated):
        return annotated, "amino_annotated"
    return plain, "plain_snp"


def build_snp_report(
    varscan_dir: str,
    sample_id: str,
    output_dir: str,
    total_reads: int = 0,
) -> str:
    os.makedirs(output_dir, exist_ok=True)
    input_path, source_label = _pick_input(varscan_dir, sample_id)
    vs_rows = _read_varscan_rows(input_path)

    out_rows: list[dict] = []
    for vs in vs_rows:
        out_rows.append(_build_row(sample_id, vs))

    out_rows.sort(key=lambda r: (r["chrom"], _try_int(r["position"]), r["alt"]))

    # Drop sites where informative depth is below MIN_TOTAL_INFORMATIVE_FRAC
    # of total mapped reads (primer-end / deletion-gap artefacts).
    if total_reads > 0:
        min_depth = MIN_TOTAL_INFORMATIVE_FRAC * total_reads
        out_rows = [r for r in out_rows if r["total_informative"] >= min_depth]


    out_path = os.path.join(output_dir, OUTPUT_FILE)
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=SNP_REPORT_COLUMNS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    print(f"[Takeda Pipeline Message] varscan_snp_report: source={source_label} "
          f"in={len(vs_rows)} wrote={len(out_rows)} -> {out_path}")
    return out_path


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Thin VarScan SNP report (no overlap with indel pipeline).",
    )
    p.add_argument("--varscan-dir", required=True,
                   help="Directory containing <sample>.snp.vcf (and optionally <sample>.snp.amino_annotated.vcf).")
    p.add_argument("--sample-id", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--total-reads", type=int, default=0,
                   help="Total mapped reads for the sample (from BAM). Sites with "
                        "total_informative < 1%% of this value are dropped.")
    return p


def main(argv: Iterable[str] | None = None) -> int:
    args = _build_parser().parse_args(list(argv) if argv is not None else None)
    build_snp_report(
        varscan_dir=args.varscan_dir,
        sample_id=args.sample_id,
        output_dir=args.output_dir,
        total_reads=args.total_reads,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
