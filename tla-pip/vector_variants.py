#!/usr/bin/env python3

"""TLA-seq -- within-vector characterization (standalone).

Within-vector characterization using coverage, VarScan small variants,
optional Pindel calls, and existing GRIDSS vector-vector breakpoints.

For each plasmid/vector contig, from a sorted BAM (mixed host+vector
reference) and the combined GRIDSS breakpoints file:

    -> plasmid mpileup (uncapped depth, -A)         [tmp]
    -> VarScan SNP / indel / consensus               (frequency-aware)
    -> per-base linear vector coverage plot          (with GRIDSS overlay)
    [-> Pindel SI/LI/D VCFs                          (optional, default off)]

Produces sample-level outputs:
    .vector_structural_breakpoints.tsv  (all vector<->vector GRIDSS records)
    .vector_characterization_summary.tsv (one row per contig)

Runs standalone (CLI below) and as a stage of ``pipeline.py``.
"""

import argparse
import csv
import os
import sys

import numpy as np

import vector_coverage_plot as vcp
from pipeline_io import get_logger

# Flat layout: tools/ helpers import by bare name (mirrors fusion_detection).
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools"))
import vector_variant_utils as vvu  # noqa: E402

log = get_logger()

_VECTOR_SV_TYPE = "vector_sv"
_NLAIII_COL = "nlaiii_category"
_STRUCT_BP_COLUMNS = (
    "contig_a", "pos_a", "contig_b", "pos_b", "breakend", "QUAL",
    "VF", "SR", "RP", "ASSR", "ASRP", "HOMLEN", "HOMSEQ",
    "n_beid", "IMPRECISE", "SB", "EVENT", _NLAIII_COL,
)
_SUMMARY_HEADER = (
    "vector_contig", "vector_length", "mean_depth", "median_depth", "max_depth",
    "n_varscan_snvs", "n_varscan_indels", "n_pindel_calls",
    "n_gridss_vector_breakpoints", "n_gridss_both_ends_nlaiii",
    "n_gridss_one_end_nlaiii", "n_gridss_neither_end_nlaiii",
)

# --- VarScan TSV parsing (adapted from ko-pip/varscan_snp_report.py) ---

_SNV_REPORT_COLUMNS = (
    "contig", "position", "ref", "alt",
    "ref_support", "alt_support", "total_informative", "total_depth",
    "var_freq", "p_value",
)
_INDEL_REPORT_COLUMNS = (
    "contig", "position", "ref", "alt", "type", "length",
    "ref_support", "alt_support", "var_freq", "p_value",
)


def _try_int(x):
    try:
        return int(x)
    except (TypeError, ValueError):
        return 0


def _parse_var_freq(token):
    """VarScan emits VarFreq as '12.34%'. Return a fraction in [0, 1]."""
    s = (token or "").strip()
    if not s:
        return 0.0
    if s.endswith("%"):
        try:
            return float(s[:-1]) / 100.0
        except ValueError:
            return 0.0
    try:
        return float(s)
    except ValueError:
        return 0.0


def _parse_varscan_tsv(path):
    """Read a VarScan pileup2snp/pileup2indel TSV; return list of dicts."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_parsed_report(rows, header, path):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(col, "")) for col in header) + "\n")
    os.replace(tmp, path)
    return len(rows)


def _read_mpileup_depths(mpileup_path):
    """Return {position_str: raw_depth} from an mpileup file."""
    depths = {}
    if not mpileup_path or not os.path.exists(mpileup_path):
        return depths
    with open(mpileup_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 4:
                depths[parts[1]] = parts[3]
    return depths


def _build_snv_report(raw_snp_path, out_path, mpileup_path=None):
    """Parse VarScan pileup2snp TSV into a clean SNV report."""
    depths = _read_mpileup_depths(mpileup_path)
    rows = []
    for vs in _parse_varscan_tsv(raw_snp_path):
        ref_support = _try_int(vs.get("Reads1", 0))
        alt_support = _try_int(vs.get("Reads2", 0))
        pos = vs.get("Position", "")
        rows.append({
            "contig": vs.get("Chrom", ""),
            "position": pos,
            "ref": vs.get("Ref", ""),
            "alt": vs.get("VarAllele") or vs.get("Var", ""),
            "ref_support": ref_support,
            "alt_support": alt_support,
            "total_informative": ref_support + alt_support,
            "total_depth": depths.get(pos, ""),
            "var_freq": f"{_parse_var_freq(vs.get('VarFreq', '')):.4f}",
            "p_value": vs.get("Pvalue", ""),
        })
    return _write_parsed_report(rows, _SNV_REPORT_COLUMNS, out_path)


def _build_indel_report(raw_indel_path, out_path):
    """Parse VarScan pileup2indel TSV into a clean indel report."""
    rows = []
    for vs in _parse_varscan_tsv(raw_indel_path):
        ref = vs.get("Ref", "")
        alt = vs.get("VarAllele") or vs.get("Var", "")
        # VarScan pileup2indel encodes deletions as "-BASES" and insertions
        # as "+BASES" in VarAllele.  Do NOT use string-length comparison
        # (ref vs alt) — the +/- prefix inflates len(alt) and inverts the call.
        indel_type = "DEL" if alt.startswith("-") else "INS"
        indel_len = len(alt) - 1 if alt and alt[0] in "+-" else len(alt)
        alt_clean = alt.lstrip("+-")
        rows.append({
            "contig": vs.get("Chrom", ""),
            "position": vs.get("Position", ""),
            "ref": ref,
            "alt": alt_clean,
            "type": indel_type,
            "length": indel_len,
            "ref_support": _try_int(vs.get("Reads1", 0)),
            "alt_support": _try_int(vs.get("Reads2", 0)),
            "var_freq": f"{_parse_var_freq(vs.get('VarFreq', '')):.4f}",
            "p_value": vs.get("Pvalue", ""),
        })
    return _write_parsed_report(rows, _INDEL_REPORT_COLUMNS, out_path)


def _read_vector_sv_breakpoints(path):
    """Return vector_sv rows from the combined gridss_breakpoints.tsv as list of dicts."""
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("breakpoint_type") == _VECTOR_SV_TYPE:
                rows.append(row)
    return rows


def _count_data_rows(path):
    """Number of data lines in a TSV (total lines minus header); 0 if absent."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return 0
    with open(path) as fh:
        return max(0, sum(1 for _ in fh) - 1)


def _count_vcf_calls(vcf_paths):
    """Count non-header lines across a list of VCF files."""
    total = 0
    for p in (vcf_paths or []):
        if os.path.exists(p):
            with open(p) as fh:
                total += sum(1 for ln in fh if not ln.startswith("#"))
    return total


def read_plasmid_contigs(plasmid_contigs_path):
    """Return the list of vector/plasmid contig names (one per line)."""
    with open(plasmid_contigs_path) as fh:
        return [line.strip() for line in fh if line.strip()]


def run(sorted_bam, mixed_fa, plasmid_contigs, out_prefix,
        gridss_breakpoints=None,
        tmp_prefix=None,
        coverage_dir=None):
    """Run within-vector characterization for every plasmid contig.

    Tmp files (mpileup, raw VarScan2 outputs, per-base depth TSV) go under
    ``<tmp_prefix>/``; user-facing outputs (SNV/indel reports, coverage plots)
    go to ``out_prefix`` and optionally ``coverage_dir``. If ``coverage_dir``
    is given, vector coverage PNGs are written there (alongside host coverage
    plots) instead of under ``out_prefix``. Returns
    ``{contig: {output: path}, "_sample_outputs": {summary/structural paths}}``.
    """
    out_dir = os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    label = os.path.basename(out_prefix)

    # Tmp layout: mpileup at <tmp_prefix>.<contig>.mpileup;
    # raw VarScan2 outputs (snp/indel/cns TSVs) under <tmp_dir>/varscan2/.
    if tmp_prefix is None:
        tmp_prefix = os.path.join(out_dir or ".", "tmp", label)
    tmp_dir = os.path.dirname(tmp_prefix)
    tmp_vs2_dir = os.path.join(tmp_dir, "varscan2")
    os.makedirs(tmp_vs2_dir, exist_ok=True)

    contigs = read_plasmid_contigs(plasmid_contigs)
    if not contigs:
        raise ValueError(f"No plasmid contigs listed in {plasmid_contigs}")
    log.info(f"Vector characterization: {len(contigs)} contig(s): {', '.join(contigs)}")

    # Read GRIDSS vector-vector breakpoints once for all contigs.
    all_sv_rows = []
    if gridss_breakpoints and os.path.exists(gridss_breakpoints):
        all_sv_rows = _read_vector_sv_breakpoints(gridss_breakpoints)
        _nlcounts = {cat: sum(1 for r in all_sv_rows if r.get(_NLAIII_COL) == cat)
                     for cat in ("both_ends_near_nlaiii", "one_end_near_nlaiii",
                                 "neither_end_near_nlaiii")}
        log.info(
            f"GRIDSS vector-vector breakpoints: {len(all_sv_rows)} total -- "
            f"both_ends_nlaiii={_nlcounts['both_ends_near_nlaiii']}, "
            f"one_end_nlaiii={_nlcounts['one_end_near_nlaiii']}, "
            f"neither_end_nlaiii={_nlcounts['neither_end_near_nlaiii']}")
    elif gridss_breakpoints:
        log.warning(f"GRIDSS breakpoints file not found: {gridss_breakpoints}")

    # Write vector structural breakpoints TSV (sample-level; all vector_sv rows).
    structural_bp_path = out_prefix + ".vector_structural_breakpoints.tsv"
    _tmp = structural_bp_path + ".tmp"
    with open(_tmp, "w") as fh:
        fh.write("\t".join(_STRUCT_BP_COLUMNS) + "\n")
        for row in all_sv_rows:
            fh.write("\t".join(str(row.get(col, "")) for col in _STRUCT_BP_COLUMNS) + "\n")
    os.replace(_tmp, structural_bp_path)
    log.info(f"Vector structural breakpoints -> {structural_bp_path}")

    fai = mixed_fa + ".fai"
    results = {}
    summary_rows = []

    for contig in contigs:
        log.info(f"---- Vector characterization: {contig} ----")
        prefix = f"{label}.{contig}"
        contig_prefix = os.path.join(out_dir, prefix) if out_dir else prefix
        tmp_contig_prefix = tmp_prefix + f".{contig}"

        # 1. mpileup (-> tmp) and VarScan SNP / indel / consensus (raw -> tmp/varscan2).
        mpileup = tmp_contig_prefix + ".mpileup"
        vvu.plasmid_mpileup(sorted_bam, mixed_fa, contig, mpileup)
        snp, indel, cns = vvu.varscan_call(
            mpileup, tmp_vs2_dir, prefix,
            min_var_freq=0.01)

        # Parse raw VarScan TSVs into clean user-facing reports.
        snv_report = contig_prefix + ".snv_report.tsv"
        indel_report = contig_prefix + ".indel_report.tsv"
        n_snvs = _build_snv_report(snp, snv_report, mpileup_path=mpileup)
        n_indels = _build_indel_report(indel, indel_report)
        log.info(f"{contig}: VarScan SNVs={n_snvs}, indels={n_indels}")

        # 2. Per-base depth (-> tmp) + characterization plot.
        length = vcp.get_contig_length(fai, contig)
        depth_tsv = tmp_contig_prefix + ".vector_depth.tsv"
        vcp.plasmid_depth(sorted_bam, contig, depth_tsv)
        depth = vcp.read_depth(depth_tsv, length)

        contig_sv_rows = [
            r for r in all_sv_rows
            if r.get("contig_a") == contig or r.get("contig_b") == contig]
        # Write vector coverage PNG to coverage_dir if provided, else vector_characterization/.
        if coverage_dir:
            os.makedirs(coverage_dir, exist_ok=True)
            char_png = os.path.join(coverage_dir, f"{label}.{contig}.vector_characterization.png")
        else:
            char_png = contig_prefix + ".vector_characterization.png"
        vcp.plot_vector_coverage(depth, length, char_png, label, contig,
                                 gridss_rows=contig_sv_rows or None)
        log.info(f"{contig}: characterization plot -> {char_png}")

        # 3. Pindel: reserved for later incorporation.
        n_pindel = "."

        summary_rows.append({
            "vector_contig": contig,
            "vector_length": length,
            "mean_depth": round(float(np.mean(depth)), 2),
            "median_depth": round(float(np.median(depth)), 2),
            "max_depth": int(np.max(depth)),
            "n_varscan_snvs": n_snvs,
            "n_varscan_indels": n_indels,
            "n_pindel_calls": n_pindel,
            "n_gridss_vector_breakpoints": len(contig_sv_rows),
            "n_gridss_both_ends_nlaiii": sum(
                1 for r in contig_sv_rows if r.get(_NLAIII_COL) == "both_ends_near_nlaiii"),
            "n_gridss_one_end_nlaiii": sum(
                1 for r in contig_sv_rows if r.get(_NLAIII_COL) == "one_end_near_nlaiii"),
            "n_gridss_neither_end_nlaiii": sum(
                1 for r in contig_sv_rows if r.get(_NLAIII_COL) == "neither_end_near_nlaiii"),
        })

        out = {
            "snv_report": snv_report,
            "indel_report": indel_report,
            "vector_characterization_png": char_png,
        }
        results[contig] = out

    # Write compact vector summary.
    summary_path = out_prefix + ".vector_characterization_summary.tsv"
    _tmp = summary_path + ".tmp"
    with open(_tmp, "w") as fh:
        fh.write("\t".join(_SUMMARY_HEADER) + "\n")
        for row in summary_rows:
            fh.write("\t".join(str(row[col]) for col in _SUMMARY_HEADER) + "\n")
    os.replace(_tmp, summary_path)
    log.info(f"Vector characterization summary -> {summary_path}")

    results["_sample_outputs"] = {
        "vector_structural_breakpoints": structural_bp_path,
        "vector_characterization_summary": summary_path,
    }
    log.info("Within-vector characterization complete.")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="TLA within-vector characterization: VarScan small variants, "
                    "per-base coverage plot (with optional GRIDSS overlay), "
                    "and optional Pindel calls.")
    parser.add_argument("-b", "--bam", required=True,
                        help="Sorted BAM aligned to the mixed (host+vector) reference")
    parser.add_argument("-r", "--mixed_fa", required=True,
                        help="Mixed reference FASTA (indexed; .fai beside it)")
    parser.add_argument("-p", "--plasmid_contigs", required=True,
                        help="plasmid_contigs.txt (vector contig names; from build_reference)")
    parser.add_argument("-o", "--out_prefix", required=True,
                        help="Output path prefix (dir is created if needed)")
    args = parser.parse_args()
    paths = run(args.bam, args.mixed_fa, args.plasmid_contigs, args.out_prefix)
    for contig, outputs in paths.items():
        if contig == "_sample_outputs":
            continue
        log.info(f"{contig}:")
        for name, path in outputs.items():
            log.info(f"  {name}: {path}")
    for name, path in paths.get("_sample_outputs", {}).items():
        log.info(f"{name}: {path}")
