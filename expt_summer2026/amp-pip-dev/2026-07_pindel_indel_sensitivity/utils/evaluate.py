#!/usr/bin/env python3
"""Score Pindel detection against the simulated truth (base -e 1 VCFs).

sim_pindel emits one base pindel2vcf VCF per SV type at -e 1 (every call, AD
supporting-read field intact). This module parses those base VCFs and applies a
supporting-read threshold in Python, so the same simulated data is scored at any
threshold without re-running Pindel:

  - Experiments A/B call evaluate() at the raw-detectability CEILING (min_support
    = config.CEILING_E = 1): "can Pindel see the event at all, ignoring the
    calling threshold".
  - Experiment C imports the building blocks (parse_condition, score_manifest)
    and sweeps min_support / min_fraction to pick an operating point, using the
    AF=0 WT controls as the false-positive baseline.

Scoring is EXACT (config.POS_TOLERANCE / SIZE_TOLERANCE = 0): a call counts as a
detection only if breakpoint AND size match the truth exactly. Same-type calls
that pass the support threshold but are not exact are false positives; near-miss
calls (within config.NEARMISS_POS_WINDOW) are written out for manual inspection.

Two allele-fraction views are reported per hit: the Pindel AD fraction
(alt / (ref + alt), which C also thresholds on) and the true mpileup fraction
(alt / raw depth at the breakpoint), which exposes how far Pindel's AD
under-counts at amplicon depth.
"""

import argparse
import csv
import functools
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config  # noqa: E402


# --- VCF parsing -------------------------------------------------------------

def parse_vcf(path):
    """Parse a (pindel2vcf) VCF into a list of call dicts."""
    calls = []
    if not os.path.exists(path):
        return calls
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 8:
                continue
            info = dict(
                kv.split("=", 1) if "=" in kv else (kv, "")
                for kv in f[7].split(";")
            )
            svtype = info.get("SVTYPE", "")
            try:
                svlen = int(info.get("SVLEN", "0"))
            except ValueError:
                svlen = 0
            if svlen == 0 and svtype in ("INS", "DEL"):
                svlen = len(f[4]) - len(f[3])
            try:
                homlen = int(info.get("HOMLEN", "0"))
            except ValueError:
                homlen = 0

            ref_reads = alt_reads = 0
            if len(f) >= 10 and ":" in f[8]:
                fd = dict(zip(f[8].split(":"), f[9].split(":")))
                ad = fd.get("AD", "")
                if "," in ad:
                    try:
                        ref_reads, alt_reads = (int(x) for x in ad.split(",")[:2])
                    except ValueError:
                        ref_reads = alt_reads = 0
            ad_depth = ref_reads + alt_reads
            calls.append({
                "pos": int(f[1]),
                "svtype": svtype,
                "svlen": svlen,
                "abs_svlen": abs(svlen),
                "homlen": homlen,
                "ref_reads": ref_reads,
                "alt_reads": alt_reads,
                "ad_frac": round(alt_reads / ad_depth, 4) if ad_depth else 0.0,
            })
    return calls


def _base_vcf(pindel_dir, svtype):
    """Base -e 1 per-type VCF paths for a given SV type (no filtered files)."""
    patterns = ["*_D.vcf"] if svtype == "DEL" else ["*_SI.vcf", "*_LI.vcf"]
    paths = []
    for pat in patterns:
        paths.extend(
            p for p in glob.glob(os.path.join(pindel_dir, pat))
            if "_Filtered" not in os.path.basename(p)
            and "." not in os.path.basename(p)[:-4].split("_")[-1]  # skip *.N.vcf
        )
    return paths


@functools.lru_cache(maxsize=None)
def parse_condition(pindel_dir):
    """Return {"DEL": [calls], "INS": [calls]} of all base calls (cached)."""
    return {
        "DEL": [c for p in _base_vcf(pindel_dir, "DEL") for c in parse_vcf(p)],
        "INS": [c for p in _base_vcf(pindel_dir, "INS") for c in parse_vcf(p)],
    }


@functools.lru_cache(maxsize=None)
def _depth_map(cond_dir):
    """1-based position -> raw depth from the condition's samtools depth file."""
    hits = glob.glob(os.path.join(cond_dir, "*.depth.txt"))
    depth = {}
    if hits:
        with open(hits[0]) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) >= 3:
                    depth[int(parts[1])] = int(parts[2])
    return depth


def mpileup_af(cond_dir, pos, alt_reads):
    """alt_reads / raw depth at pos (true allele fraction). '' if unavailable."""
    depth = _depth_map(cond_dir).get(pos)
    if not depth:
        return ""
    return round(alt_reads / depth, 4)


# --- matching ----------------------------------------------------------------

def size_matches(reported, expected):
    tol = max(config.SIZE_TOLERANCE, expected * config.SIZE_TOLERANCE_FRAC)
    return abs(reported - expected) <= tol


def is_exact(call, expected_svtype, expected_len, truth_pos):
    return (call["svtype"] == expected_svtype
            and abs(call["pos"] - truth_pos) <= config.POS_TOLERANCE
            and size_matches(call["abs_svlen"], expected_len))


def passes_support(call, min_support, min_fraction):
    if min_support is not None and call["alt_reads"] < min_support:
        return False
    if min_fraction is not None and call["ad_frac"] < min_fraction:
        return False
    return True


# --- per-condition scoring ---------------------------------------------------

def score_mutant(pindel_dir, cond_dir, expected_svtype, expected_len, truth_pos,
                 min_support, min_fraction):
    """Score one mutant condition at a single threshold.

    Returns (row, nearmiss_calls). row has exact/offtarget counts and the best
    hit's reported fields (exact hit if any, else nearest near-miss diagnostic).
    """
    calls = [c for c in parse_condition(pindel_dir)[expected_svtype]
             if passes_support(c, min_support, min_fraction)]
    exact = [c for c in calls
             if is_exact(c, expected_svtype, expected_len, truth_pos)]
    offtarget = [c for c in calls if c not in exact]
    near = [c for c in offtarget
            if abs(c["pos"] - truth_pos) <= config.NEARMISS_POS_WINDOW]

    if exact:
        best = exact[0]
    elif near:
        best = min(near, key=lambda c: (abs(c["pos"] - truth_pos),
                                        abs(c["abs_svlen"] - expected_len)))
    else:
        best = None

    row = {
        "detected": int(bool(exact)),
        "n_offtarget": len(offtarget),
        "reported_pos": best["pos"] if best else "",
        "reported_svlen": best["abs_svlen"] if best else "",
        "reported_homlen": best["homlen"] if best else "",
        "pos_offset": (best["pos"] - truth_pos) if best else "",
        "alt_reads": best["alt_reads"] if best else "",
        "ref_reads": best["ref_reads"] if best else "",
        "ad_frac": best["ad_frac"] if best else "",
        "mpileup_af": mpileup_af(cond_dir, truth_pos, best["alt_reads"]) if best else "",
    }
    return row, near


def score_control(pindel_dir, min_support, min_fraction):
    """Spurious DEL/INS call counts in a WT-only control at one threshold."""
    cond = parse_condition(pindel_dir)
    return {
        "DEL": sum(1 for c in cond["DEL"]
                   if passes_support(c, min_support, min_fraction)),
        "INS": sum(1 for c in cond["INS"]
                   if passes_support(c, min_support, min_fraction)),
    }


# --- run table ---------------------------------------------------------------

def load_manifest(outdir, experiment):
    """Read one experiment's run table (results/manifest.csv)."""
    manifest_path = os.path.join(outdir, "results", "manifest.csv")
    if not os.path.exists(manifest_path):
        return None
    with open(manifest_path) as fh:
        return list(csv.DictReader(fh))


def _is_control(m):
    return str(m.get("is_control", "")).strip() in ("1", "True", "true")


def score_manifest(outdir, experiment, min_support=config.CEILING_E,
                   min_fraction=None):
    """Score every condition in an experiment's run table at one threshold.

    Returns (results_df, mutant_records, control_fp, nearmiss_rows):
      results_df      per-condition detail (mutants + controls)
      mutant_records  flat rows (group keys + exact/offtarget) for confusion math
      control_fp      {(dataset, read_len, depth): {"DEL": n, "INS": n}}
      nearmiss_rows   same-type near-miss calls for manual inspection
    """
    manifest = load_manifest(outdir, experiment)
    if manifest is None:
        raise FileNotFoundError(
            f"No manifest at {outdir}/results/manifest.csv - run "
            f"the experiment first.")

    rows, mutant_records, nearmiss_rows = [], [], []
    control_fp = {}
    for m in manifest:
        cond_id = m["condition_id"]
        cond_dir = os.path.join(outdir, "conditions", cond_id)
        pindel_dir = os.path.join(cond_dir, "pindel-output")
        read_len = int(m["read_len"])
        depth = int(m["depth"])
        seed = int(m["seed"])
        dataset = m.get("dataset", "")
        experiment = m.get("experiment", "")

        base = {
            "condition_id": cond_id, "experiment": experiment,
            "dataset": dataset, "is_control": 1 if _is_control(m) else 0,
            "read_len": read_len, "depth": depth, "seed": seed,
        }

        if _is_control(m):
            fp = score_control(pindel_dir, min_support, min_fraction)
            control_fp.setdefault((dataset, read_len, depth), {}).update(
                {"DEL": fp["DEL"], "INS": fp["INS"]})
            base.update({
                "indel_type": "", "length": "", "allele_freq": 0.0,
                "amplicon_pos": "", "expected_svtype": "", "truth_pos": "",
                "detected": "", "n_offtarget": fp["DEL"] + fp["INS"],
                "fp_del": fp["DEL"], "fp_ins": fp["INS"],
            })
            rows.append(base)
            continue

        expected_svtype = m["expected_svtype"]
        length = int(m["length"])
        truth_pos = int(m["anchor_pos"])
        amplicon_pos = int(m["amplicon_pos"])
        af = float(m["allele_freq"])
        row, near = score_mutant(pindel_dir, cond_dir, expected_svtype, length,
                                 truth_pos, min_support, min_fraction)
        base.update({
            "indel_type": m["indel_type"], "length": length, "allele_freq": af,
            "amplicon_pos": amplicon_pos, "expected_svtype": expected_svtype,
            "truth_pos": truth_pos,
        })
        base.update(row)
        rows.append(base)

        mutant_records.append({
            "dataset": dataset, "indel_type": m["indel_type"], "length": length,
            "allele_freq": af, "read_len": read_len, "depth": depth, "seed": seed,
            "exact": row["detected"], "offtarget": row["n_offtarget"],
        })
        for c in near:
            nearmiss_rows.append({
                "condition_id": cond_id, "experiment": experiment,
                "dataset": dataset, "indel_type": m["indel_type"],
                "expected_svtype": expected_svtype, "expected_svlen": length,
                "truth_pos": truth_pos, "read_len": read_len, "depth": depth,
                "seed": seed, "reported_pos": c["pos"],
                "reported_svlen": c["abs_svlen"],
                "pos_offset": c["pos"] - truth_pos,
                "size_offset": c["abs_svlen"] - length,
                "homlen": c["homlen"], "alt_reads": c["alt_reads"],
                "ad_frac": c["ad_frac"],
                "mpileup_af": mpileup_af(cond_dir, truth_pos, c["alt_reads"]),
            })

    return pd.DataFrame(rows), mutant_records, control_fp, nearmiss_rows


# --- aggregation -------------------------------------------------------------

def _safe_div(num, den):
    return num / den if den else ""


def confusion_metrics(mutant_records, control_fp):
    """Per-group TP/FP/FN + sensitivity/precision/F1 pooled over seeds."""
    df = pd.DataFrame(mutant_records)
    if df.empty:
        return pd.DataFrame()
    group_cols = ["dataset", "indel_type", "length", "allele_freq",
                  "read_len", "depth"]
    rows = []
    for keys, g in df.groupby(group_cols):
        dataset, indel_type, length, af, read_len, depth = keys
        svtype = "DEL" if indel_type == "del" else "INS"
        tp = int(g["exact"].sum())
        n = len(g)
        fn = n - tp
        fp = int(g["offtarget"].sum())
        fp += int(control_fp.get((dataset, read_len, depth), {}).get(svtype, 0))
        sens = _safe_div(tp, tp + fn)
        prec = _safe_div(tp, tp + fp)
        f1 = ""
        if prec != "" and sens != "" and (prec + sens) > 0:
            f1 = round(2 * prec * sens / (prec + sens), 4)
        rows.append({
            "dataset": dataset, "indel_type": indel_type, "length": length,
            "allele_freq": af, "read_len": read_len, "depth": depth,
            "n": n, "TP": tp, "FP": fp, "FN": fn,
            "sensitivity": round(sens, 4) if sens != "" else "",
            "precision": round(prec, 4) if prec != "" else "",
            "f1": f1,
        })
    return pd.DataFrame(rows).sort_values(group_cols)


def _pivot_sensitivity(results, out_path):
    """Exp A exact-detection pivot: rows=length, cols=AF, per dataset/type/rl/depth."""
    sens = (results.groupby(["dataset", "indel_type", "read_len", "depth",
                             "length", "allele_freq"])["detected"]
            .mean().reset_index())
    with open(out_path, "w") as fh:
        fh.write("# Pindel EXACT detection rate at the -e ceiling, "
                 "averaged over seeds\n")
        for dataset in sorted(sens["dataset"].unique()):
            for indel_type in sorted(sens["indel_type"].unique()):
                for read_len in sorted(sens["read_len"].unique()):
                    for depth in sorted(sens["depth"].unique()):
                        sub = sens[(sens["dataset"] == dataset)
                                   & (sens["indel_type"] == indel_type)
                                   & (sens["read_len"] == read_len)
                                   & (sens["depth"] == depth)]
                        if sub.empty:
                            continue
                        pivot = sub.pivot(index="length", columns="allele_freq",
                                          values="detected")
                        fh.write(f"# dataset={dataset} indel_type={indel_type} "
                                 f"read_len={read_len} depth={depth} "
                                 f"(rows=length bp, cols=allele freq)\n")
                        pivot.to_csv(fh)
                        fh.write("\n")


def _pivot_position(results, out_path):
    """Exp B exact-detection pivot: rows=breakpoint position, cols=type."""
    sens = (results.groupby(["dataset", "read_len", "indel_type",
                             "amplicon_pos"])["detected"]
            .mean().reset_index())
    with open(out_path, "w") as fh:
        fh.write("# Exp B exact-detection rate vs breakpoint position "
                 "at the -e ceiling, averaged over seeds\n")
        for dataset in sorted(sens["dataset"].unique()):
            for read_len in sorted(sens["read_len"].unique()):
                sub = sens[(sens["dataset"] == dataset)
                           & (sens["read_len"] == read_len)]
                if sub.empty:
                    continue
                pivot = sub.pivot(index="amplicon_pos", columns="indel_type",
                                  values="detected")
                fh.write(f"# dataset={dataset} read_len={read_len} "
                         f"(rows=amplicon position, cols=indel type)\n")
                pivot.to_csv(fh)
                fh.write("\n")


def evaluate(outdir, experiment, min_support=config.CEILING_E,
             min_fraction=None):
    """Score one experiment at a single threshold and write its CSV suite.

    Writes results/ detection_results.csv, and (Exp A)
    sensitivity_by_length_af.csv, confusion_metrics.csv, nearmiss.csv, or (Exp B)
    position_sensitivity.csv. Returns the per-condition results DataFrame.
    """
    results, mutant_records, control_fp, nearmiss_rows = score_manifest(
        outdir, experiment, min_support, min_fraction)
    results_dir = os.path.join(outdir, "results")
    os.makedirs(results_dir, exist_ok=True)

    results.to_csv(os.path.join(results_dir, "detection_results.csv"),
                   index=False)

    exp_a = results[(results["experiment"] == "A") & (results["is_control"] == 0)]
    if not exp_a.empty:
        _pivot_sensitivity(
            exp_a, os.path.join(results_dir, "sensitivity_by_length_af.csv"))

    conf = confusion_metrics(mutant_records, control_fp)
    if not conf.empty:
        conf.to_csv(os.path.join(results_dir, "confusion_metrics.csv"),
                    index=False)

    if nearmiss_rows:
        pd.DataFrame(nearmiss_rows).to_csv(
            os.path.join(results_dir, "nearmiss.csv"), index=False)

    exp_b = results[(results["experiment"] == "B") & (results["is_control"] == 0)]
    if not exp_b.empty:
        _pivot_position(
            exp_b, os.path.join(results_dir, "position_sensitivity.csv"))

    return results


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default=config.OUTDIR)
    parser.add_argument("--experiment", required=True,
                        help="Experiment tag whose manifest to score (A / B)")
    parser.add_argument("--min-support", type=int, default=config.CEILING_E,
                        help="Absolute alt-supporting-read floor (default: "
                             "%(default)s = raw-detectability ceiling)")
    args = parser.parse_args()

    results = evaluate(args.outdir, args.experiment,
                       min_support=args.min_support)
    mutants = results[results["is_control"] == 0] if not results.empty else results
    n = len(mutants)
    detected = int(pd.to_numeric(mutants["detected"], errors="coerce")
                   .fillna(0).sum()) if n else 0
    print(f"[evaluate] EXACT detections: {detected}/{n} mutant conditions "
          f"(min_support={args.min_support})")


if __name__ == "__main__":
    main()
