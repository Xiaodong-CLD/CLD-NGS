#!/usr/bin/env python3
"""Experiment C - choosing the pindel2vcf support threshold (-e).

QUESTIONS ADDRESSED
  1. Given the raw calls Pindel produces, what supporting-read threshold best
     separates true low-fraction indels from false positives? i.e. where should
     pindel2vcf's -e be set for amplicon data?
  2. Is an ABSOLUTE read count (-e N) or a FRACTIONAL cutoff (alt / (ref+alt))
     the better operating rule across depths?
  3. How badly does Pindel's AD supporting-read count under-report the true
     allele fraction at amplicon depth (Pindel AD fraction vs the mpileup
     fraction from raw depth)? This says whether an AD-based threshold is even
     trustworthy here.

This is a POST-HOC re-score: no new reads are simulated. It reuses the base
-e 1 VCFs already produced by Experiment A (run detect_sensitivity_grid.py
first) and sweeps the threshold in Python over a challenging subset (low allele
fraction, short-to-moderate indels) plus Experiment A's WT controls for the
false-positive baseline.
"""

import argparse
import os
import sys

import pandas as pd

BENCH_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH_DIR)
sys.path.insert(0, os.path.join(BENCH_DIR, "utils"))

import config     # noqa: E402
import evaluate   # noqa: E402

SOURCE_EXPERIMENT = "A"   # whose base VCFs / manifest we re-score
EXPERIMENT = "C"

# --- Experiment C sweep (owned here, not in config) --------------------------
# Absolute supporting-read floors (pindel2vcf -e) and fractional cutoffs.
# The absolute ladder is extended past 100 so the F1 peak stays bracketed at the
# high-depth tiers (the -e ~= 0.001 x depth rule predicts optima near 5 / 10 /
# 100 for 5000 / 10000 / 100000 depth).
E_ABSOLUTE = [1, 3, 5, 10, 15, 20, 30, 50, 75, 100, 150, 200, 300, 500]
E_FRACTION = [0.005, 0.01, 0.02, 0.05]
# Challenging subset of Experiment A: low allele fraction, short/moderate indels.
# Every depth tier present in Experiment A is swept independently so the optimal
# threshold can be compared across depths.
SUBSET_AF = [0.01, 0.02]
SUBSET_LENGTHS = [1, 5, 10, 20]
SUBSET_DEPTHS = [5000, 10000, 100000]


def _cond_paths(source_outdir, cond_id):
    cond_dir = os.path.join(source_outdir, "conditions", cond_id)
    return cond_dir, os.path.join(cond_dir, "pindel-output")


def _select(manifest, depth):
    """Split the source manifest into the challenging mutant subset + controls
    for a single depth tier."""
    mutants, controls = [], []
    for m in manifest:
        if int(m["depth"]) != depth:
            continue
        if evaluate._is_control(m):
            controls.append(m)
        elif (float(m["allele_freq"]) in SUBSET_AF
              and int(m["length"]) in SUBSET_LENGTHS):
            mutants.append(m)
    return mutants, controls


def _score_at(source_outdir, mutants, controls, min_support, min_fraction):
    """Per-indel-type TP/FN/FP over the subset at one threshold.

    Returns dict mapping indel_type -> (tp, fn, fp), including a pooled 'all'.
    """
    counts = {}  # indel_type -> [tp, fn, fp]
    for m in mutants:
        cond_dir, pindel_dir = _cond_paths(source_outdir, m["condition_id"])
        row, _ = evaluate.score_mutant(
            pindel_dir, cond_dir, m["expected_svtype"], int(m["length"]),
            int(m["anchor_pos"]), min_support, min_fraction)
        itype = m["indel_type"]
        if itype not in counts:
            counts[itype] = [0, 0, 0]
        counts[itype][0] += row["detected"]
        counts[itype][1] += 1 - row["detected"]
        counts[itype][2] += row["n_offtarget"]
    for c in controls:
        _, pindel_dir = _cond_paths(source_outdir, c["condition_id"])
        spurious = evaluate.score_control(pindel_dir, min_support, min_fraction)
        for itype, key in [("del", "DEL"), ("ins", "INS")]:
            if itype not in counts:
                counts[itype] = [0, 0, 0]
            counts[itype][2] += spurious[key]
    # Pooled across types
    counts["all"] = [
        sum(v[0] for v in counts.values()),
        sum(v[1] for v in counts.values()),
        sum(v[2] for v in counts.values()),
    ]
    return counts


def _metrics(depth, threshold_type, value, tp, fn, fp, indel_type="all"):
    sens = tp / (tp + fn) if (tp + fn) else ""
    prec = tp / (tp + fp) if (tp + fp) else ""
    f1 = ""
    if sens != "" and prec != "" and (sens + prec) > 0:
        f1 = round(2 * sens * prec / (sens + prec), 4)
    return {
        "depth": depth,
        "indel_type": indel_type,
        "threshold_type": threshold_type, "threshold": value,
        "n_mutant": tp + fn, "TP": tp, "FN": fn, "FP": fp,
        "sensitivity": round(sens, 4) if sens != "" else "",
        "precision": round(prec, 4) if prec != "" else "",
        "f1": f1,
    }


def _af_gap(source_outdir, mutants):
    """Per-condition Pindel AD fraction vs true mpileup fraction at the ceiling."""
    rows = []
    for m in mutants:
        cond_dir, pindel_dir = _cond_paths(source_outdir, m["condition_id"])
        row, _ = evaluate.score_mutant(
            pindel_dir, cond_dir, m["expected_svtype"], int(m["length"]),
            int(m["anchor_pos"]), config.CEILING_E, None)
        if row["detected"]:
            rows.append({
                "condition_id": m["condition_id"], "dataset": m["dataset"],
                "indel_type": m["indel_type"], "length": int(m["length"]),
                "allele_freq": float(m["allele_freq"]),
                "read_len": int(m["read_len"]), "depth": int(m["depth"]),
                "seed": int(m["seed"]), "alt_reads": row["alt_reads"],
                "ref_reads": row["ref_reads"], "ad_frac": row["ad_frac"],
                "mpileup_af": row["mpileup_af"],
            })
    return rows


def _score_dataset(dataset):
    """Re-score Experiment A's base VCFs for one dataset; write Exp C CSVs.

    Sweeps every depth tier in SUBSET_DEPTHS that is present in the manifest.
    Returns True if the dataset was scored, False if it was skipped (no
    Experiment A manifest, or no depth tier had a non-empty subset).
    """
    source_outdir = os.path.join(config.OUTDIR, f"expt{SOURCE_EXPERIMENT}", dataset)
    manifest = evaluate.load_manifest(source_outdir, SOURCE_EXPERIMENT)
    if manifest is None:
        print(f"[C] {dataset}: no Experiment {SOURCE_EXPERIMENT} manifest at "
              f"{source_outdir}/results/manifest.csv - run "
              f"detect_sensitivity_grid.py for this dataset first; skipping.")
        return False

    sweep, gap, tiers = [], [], []
    for depth in SUBSET_DEPTHS:
        mutants, controls = _select(manifest, depth)
        if not mutants:
            continue
        for e in E_ABSOLUTE:
            counts = _score_at(source_outdir, mutants, controls, e, None)
            for itype, (tp, fn, fp) in counts.items():
                sweep.append(_metrics(depth, "absolute", e, tp, fn, fp, itype))
        for frac in E_FRACTION:
            counts = _score_at(source_outdir, mutants, controls, None, frac)
            for itype, (tp, fn, fp) in counts.items():
                sweep.append(_metrics(depth, "fraction", frac, tp, fn, fp, itype))
        gap.extend(_af_gap(source_outdir, mutants))
        tiers.append((depth, len(mutants), len(controls)))

    if not tiers:
        print(f"[C] {dataset}: challenging subset is empty at every depth - "
              f"check SUBSET_AF / SUBSET_LENGTHS / SUBSET_DEPTHS against the "
              f"Experiment A grid; skipping.")
        return False

    results_dir = os.path.join(config.OUTDIR, f"expt{EXPERIMENT}", dataset, "results")
    os.makedirs(results_dir, exist_ok=True)
    pd.DataFrame(sweep).to_csv(
        os.path.join(results_dir, "threshold_sweep.csv"), index=False)
    pd.DataFrame(gap).to_csv(
        os.path.join(results_dir, "af_gap.csv"), index=False)

    tier_str = ", ".join(f"d{d}: {nm}mut/{nc}ctrl" for d, nm, nc in tiers)
    print(f"[C] {dataset}: swept {len(E_ABSOLUTE)} absolute + {len(E_FRACTION)} "
          f"fractional thresholds over {len(tiers)} depth tier(s) [{tier_str}] "
          f"-> {results_dir}")
    return True


def main(datasets=None):
    if datasets is None:
        datasets = list(config.DATASETS)
    scored = [d for d in datasets if _score_dataset(d)]
    if not scored:
        sys.exit(f"[C] No datasets scored - run detect_sensitivity_grid.py "
                 f"(Experiment {SOURCE_EXPERIMENT}) first.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Experiment C: support-threshold sweep over Experiment A VCFs.")
    parser.add_argument("-d", "--dataset", choices=list(config.DATASETS),
                        help="Score one dataset only (omit to score all).")
    args = parser.parse_args()
    main(datasets=[args.dataset] if args.dataset else None)
