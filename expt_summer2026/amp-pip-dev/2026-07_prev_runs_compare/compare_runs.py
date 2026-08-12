#!/usr/bin/env python3
"""Compare amp-pip-dev VarScan variant calls between two pipeline runs.

For every sample present in *both* Results directories, compare SNP and indel
VarScan VCFs position-by-position.  Reports per-sample:
  - MATCHED   : same (Position, Ref, VarAllele) in both runs + VarFreq delta
  - OLD-ONLY  : variant called in old but not new
  - NEW-ONLY  : variant called in new but not old

Writes per-sample TSVs and a summary to --output-dir.

Usage (from the 2026-07_prev_runs_compare/ directory):
  python compare_runs.py \\
      --old-results /mnt/z/.../Results-1 \\
      --new-results AvB1_HC/Results/ \\
      --output-dir comparison/AvB1_HC/

Reconstructed Aug 2026 from shell history after accidental deletion.
"""
import argparse
import csv
import os
import sys


# ── VarScan parser ────────────────────────────────────────────────────────────

def parse_varscan(path):
    """Return list of dicts from a VarScan pileup2snp / pileup2indel file.

    Each dict has at least: Position(int), Ref, VarAllele, VarFreq(str with %).
    Returns empty list if file missing.
    """
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                continue
            row = dict(zip(header, fields))
            try:
                row["Position"] = int(row["Position"])
            except (KeyError, ValueError):
                continue
            rows.append(row)
    return rows


def _freq(s):
    """'4.23%' -> 4.23"""
    try:
        return float(s.rstrip("%"))
    except (ValueError, AttributeError, TypeError):
        return None


def _key(row):
    return (row["Position"], row.get("Ref", ""), row.get("VarAllele", row.get("Var", "")))


# ── Comparison logic ──────────────────────────────────────────────────────────

def compare_calls(old_rows, new_rows):
    """Return (matched, old_only, new_only) lists."""
    old_map = {_key(r): r for r in old_rows}
    new_map = {_key(r): r for r in new_rows}
    okeys = set(old_map)
    nkeys = set(new_map)

    matched = []
    for k in sorted(okeys & nkeys):
        old_freq = _freq(old_map[k].get("VarFreq", ""))
        new_freq = _freq(new_map[k].get("VarFreq", ""))
        delta = None
        if old_freq is not None and new_freq is not None:
            delta = new_freq - old_freq
        matched.append({
            "Position": k[0], "Ref": k[1], "VarAllele": k[2],
            "OldFreq": old_map[k].get("VarFreq", ""),
            "NewFreq": new_map[k].get("VarFreq", ""),
            "DeltaFreq": f"{delta:+.2f}pp" if delta is not None else ".",
        })

    old_only = [
        {"Position": k[0], "Ref": k[1], "VarAllele": k[2],
         "VarFreq": old_map[k].get("VarFreq", "")}
        for k in sorted(okeys - nkeys)
    ]
    new_only = [
        {"Position": k[0], "Ref": k[1], "VarAllele": k[2],
         "VarFreq": new_map[k].get("VarFreq", "")}
        for k in sorted(nkeys - okeys)
    ]
    return matched, old_only, new_only


def _write_tsv(path, rows, fieldnames):
    if not rows:
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# ── Per-sample comparison ─────────────────────────────────────────────────────

def sample_dirs(results_dir):
    """Return sorted list of sample names (immediate subdirectories)."""
    if not os.path.isdir(results_dir):
        return []
    return sorted(
        d for d in os.listdir(results_dir)
        if os.path.isdir(os.path.join(results_dir, d))
    )


def compare_sample(sample, old_results, new_results, output_dir):
    """Compare VarScan SNP + indel calls for one sample. Returns summary dict."""
    summary = {"sample": sample}
    for kind in ("snp", "indel"):
        old_path = os.path.join(old_results, sample, "varscan-output",
                                f"{sample}.{kind}.vcf")
        new_path = os.path.join(new_results, sample, "varscan-output",
                                f"{sample}.{kind}.vcf")
        old_rows = parse_varscan(old_path)
        new_rows = parse_varscan(new_path)
        matched, old_only, new_only = compare_calls(old_rows, new_rows)

        summary[f"{kind}_old"] = len(old_rows)
        summary[f"{kind}_new"] = len(new_rows)
        summary[f"{kind}_matched"] = len(matched)
        summary[f"{kind}_old_only"] = len(old_only)
        summary[f"{kind}_new_only"] = len(new_only)

        prefix = os.path.join(output_dir, f"{sample}.{kind}")
        if matched:
            _write_tsv(f"{prefix}.matched.tsv", matched,
                       ["Position", "Ref", "VarAllele", "OldFreq", "NewFreq", "DeltaFreq"])
        if old_only:
            _write_tsv(f"{prefix}.old_only.tsv", old_only,
                       ["Position", "Ref", "VarAllele", "VarFreq"])
        if new_only:
            _write_tsv(f"{prefix}.new_only.tsv", new_only,
                       ["Position", "Ref", "VarAllele", "VarFreq"])

        print(f"  [{kind}] old={len(old_rows)} new={len(new_rows)} "
              f"matched={len(matched)} old_only={len(old_only)} new_only={len(new_only)}")

    return summary


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Compare VarScan calls between two amp-pip-dev runs.")
    parser.add_argument("--old-results", required=True,
                        help="Path to the old (baseline) Results/ directory.")
    parser.add_argument("--new-results", required=True,
                        help="Path to the new Results/ directory.")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for comparison TSVs and summary.")
    args = parser.parse_args()

    old_results = args.old_results
    new_results = args.new_results
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    old_samples = set(sample_dirs(old_results))
    new_samples = set(sample_dirs(new_results))
    common = sorted(old_samples & new_samples)

    if not common:
        print(f"No common samples between {old_results} and {new_results}.",
              file=sys.stderr)
        return 1

    print(f"Comparing {len(common)} common sample(s)  "
          f"(old_only={len(old_samples - new_samples)}  "
          f"new_only={len(new_samples - old_samples)})")

    summaries = []
    for sample in common:
        print(f"\n=== {sample} ===")
        summaries.append(
            compare_sample(sample, old_results, new_results, output_dir))

    # Write summary TSV
    summary_path = os.path.join(output_dir, "summary.tsv")
    fields = ["sample",
              "snp_old", "snp_new", "snp_matched", "snp_old_only", "snp_new_only",
              "indel_old", "indel_new", "indel_matched", "indel_old_only", "indel_new_only"]
    _write_tsv(summary_path, summaries, fields)
    print(f"\nSummary written to {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
