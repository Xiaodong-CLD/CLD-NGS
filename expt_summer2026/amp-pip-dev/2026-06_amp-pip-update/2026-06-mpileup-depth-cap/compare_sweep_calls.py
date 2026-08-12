#!/usr/bin/env python3
"""compare_sweep_calls.py
Compare variant CALLS from the mpileup -d sweep against d=0 (truth).

For each finite -d value and each sample, classifies called positions
(VarFreq >= CALL_THRESHOLD) as:
  shared   – called at both d=0 and finite d
  old_only – called at d=0 but dropped below threshold at finite d
  new_only – called at finite d but absent from d=0

Breaks down by frequency bin (using d=0 freq for shared/old_only; finite-d
freq for new_only):
  ge2   : VarFreq >= 2%
  1to2  : 1% <= VarFreq < 2%

Writes:
  <output_dir>/sweep_variant_comparison.tsv  – one row per variant × (sample, d)
  <output_dir>/sweep_call_summary.tsv        – counts per (sample, d)

Usage:
    python compare_sweep_calls.py <run_dir> [--output-dir <dir>]
"""

import argparse
import glob
import os
import sys

import pandas as pd

D_VALUES = [0, 8000, 12000, 16000]
REF_D = 0
CALL_THRESHOLD = 1.0   # percent; positions with VarFreq >= this are "called"

CONSENSUS_COLS = [
    "Chrom", "Position", "Ref", "Cons", "Reads1", "Reads2", "VarFreq",
    "Strands1", "Strands2", "Qual1", "Qual2", "Pvalue",
    "MapQual1", "MapQual2", "Reads1Plus", "Reads1Minus",
    "Reads2Plus", "Reads2Minus", "VarAllele",
]
MATCH_KEY = ["Chrom", "Position", "Ref", "VarAllele"]


def parse_consensus(path, sample):
    """Read a VarScan pileup2cns consensus file.

    Returns a DataFrame of called positions (VarFreq >= CALL_THRESHOLD)
    with columns: SampleID, Chrom, Position, Ref, VarAllele, VarFreq_pct.
    """
    df = pd.read_csv(path, sep="\t", header=0, names=CONSENSUS_COLS,
                     usecols=["Chrom", "Position", "Ref", "VarAllele", "VarFreq"])
    df["VarFreq_pct"] = df["VarFreq"].str.rstrip("%").astype(float)
    df = df[df["VarFreq_pct"] >= CALL_THRESHOLD].copy()
    df["SampleID"] = sample
    return df[["SampleID", "Chrom", "Position", "Ref", "VarAllele", "VarFreq_pct"]]


def discover_samples(run_dir):
    """Find samples with all expected sweep consensus files."""
    samples = []
    for sweep_dir in sorted(glob.glob(os.path.join(run_dir, "Results", "*", "depth-sweep"))):
        sample = os.path.basename(os.path.dirname(sweep_dir))
        paths = {d: os.path.join(sweep_dir, f"{sample}.d{d}.consensus.vcf")
                 for d in D_VALUES}
        if all(os.path.isfile(p) for p in paths.values()):
            samples.append((sample, paths))
    return samples


def freq_bin(vf):
    """Return 'ge2' for VarFreq >= 2%, '1to2' for 1-2%."""
    return "ge2" if vf >= 2.0 else "1to2"


def compare_one(ref_df, new_df, d, sample):
    """Merge reference (d=0) vs finite-d calls; return annotated comparison DataFrame."""
    key = ["SampleID"] + MATCH_KEY
    ref = ref_df[key + ["VarFreq_pct"]].rename(columns={"VarFreq_pct": "ref_freq"})
    new = new_df[key + ["VarFreq_pct"]].rename(columns={"VarFreq_pct": "new_freq"})

    merged = pd.merge(ref, new, on=key, how="outer", indicator=True)
    merged["status"] = merged["_merge"].map(
        {"both": "shared", "left_only": "old_only", "right_only": "new_only"}
    )
    # Frequency used for binning: reference freq for shared/old_only, finite-d for new_only
    merged["bin_freq"] = merged.apply(
        lambda r: r["new_freq"] if r["status"] == "new_only" else r["ref_freq"], axis=1
    )
    merged["freq_bin"] = merged["bin_freq"].apply(freq_bin)
    merged["d"] = d
    merged["sample"] = sample
    return merged.drop(columns=["_merge", "bin_freq"])


def build_summary(cmp_df):
    """Pivot comparison DataFrame into a counts table per (sample, d)."""
    rows = []
    for (sample, d), grp in cmp_df.groupby(["sample", "d"]):
        row = {"sample": sample, "d": d}
        for status in ("shared", "old_only", "new_only"):
            sub = grp[grp["status"] == status]
            row[f"{status}_ge2"]  = (sub["freq_bin"] == "ge2").sum()
            row[f"{status}_1to2"] = (sub["freq_bin"] == "1to2").sum()
        # Totals for convenience
        row["old_total"] = row["shared_ge2"] + row["shared_1to2"] + row["old_only_ge2"] + row["old_only_1to2"]
        row["new_total"] = row["shared_ge2"] + row["shared_1to2"] + row["new_only_ge2"] + row["new_only_1to2"]
        rows.append(row)
    cols = [
        "sample", "d",
        "old_total", "new_total",
        "shared_ge2", "shared_1to2",
        "old_only_ge2", "old_only_1to2",
        "new_only_ge2", "new_only_1to2",
    ]
    return pd.DataFrame(rows, columns=cols)


def main():
    parser = argparse.ArgumentParser(description="Compare sweep variant calls vs d=0.")
    parser.add_argument("run_dir", help="Root run directory (contains Results/)")
    parser.add_argument("--output-dir", default=None,
                        help="Directory for output TSVs (default: <run_dir>/depth-sweep-comparison)")
    args = parser.parse_args()

    if not os.path.isdir(args.run_dir):
        sys.exit(f"ERROR: run_dir not found: {args.run_dir}")

    output_dir = args.output_dir or os.path.join(args.run_dir, "depth-sweep-comparison")
    os.makedirs(output_dir, exist_ok=True)

    samples = discover_samples(args.run_dir)
    if not samples:
        sys.exit(f"ERROR: no sweep consensus files found under {args.run_dir}/Results/*/depth-sweep/")

    print(f"Found {len(samples)} sample(s): {[s for s, _ in samples]}")
    print(f"Call threshold: {CALL_THRESHOLD}%  |  d=0 is reference (truth)\n")

    all_comparisons = []
    for sample, paths in samples:
        ref_df = parse_consensus(paths[REF_D], sample)
        print(f"  {sample}: d=0 → {len(ref_df)} called variants")
        for d in [v for v in D_VALUES if v != REF_D]:
            new_df = parse_consensus(paths[d], sample)
            cmp = compare_one(ref_df, new_df, d, sample)
            all_comparisons.append(cmp)

    all_cmp = pd.concat(all_comparisons, ignore_index=True)
    summary = build_summary(all_cmp)

    cmp_path = os.path.join(output_dir, "sweep_variant_comparison.tsv")
    sum_path = os.path.join(output_dir, "sweep_call_summary.tsv")
    all_cmp.drop(columns=["ref_freq", "new_freq"], errors="ignore").to_csv(cmp_path, sep="\t", index=False)
    summary.to_csv(sum_path, sep="\t", index=False)

    print(f"\nWrote: {cmp_path}")
    print(f"Wrote: {sum_path}")
    print("\nSummary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
