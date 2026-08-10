#!/usr/bin/env python3
"""Experiment D: head-to-head indel-caller sensitivity comparison.

Re-uses Experiment A's already-simulated condition BAMs (no re-simulation) and
re-calls indels with several callers, scoring each against the SAME simulated
truth with the SAME exact matching used by Experiment A (evaluate.is_exact,
config.POS_TOLERANCE / SIZE_TOLERANCE = 0). This isolates caller sensitivity
from the simulation, so Pindel / VarScan2 / Mutect2 / LoFreq are compared on
identical input.

Callers
-------
  pindel   : incumbent. Re-uses Exp A's base -e 1 pindel-output VCFs, but
             re-scored at a fractional cutoff on Pindel's own AD
             (alt / (alt + ref) >= PINDEL_OPT_FRACTION), the operating rule the
             Exp C sweep found transferable across depth. Numbers here reflect
             that operating point rather than Exp A's -e 1 ceiling. Note the
             cutoff is a fraction of Pindel AD (which overstates true VAF ~4x,
             see Exp C), not a fraction of true depth.
  varscan  : samtools mpileup | varscan mpileup2indel (pileup caller).
  mutect2  : GATK Mutect2 tumor-only (local reassembly, low-AF/somatic).
  lofreq   : LoFreq call --call-indels (quality-aware, low-AF).

Why these settings
------------------
  * Depth is held at a small fixed tier (DEPTHS below); the caller comparison
    asks "who detects the event", not "how does detection scale with depth"
    (that is Exp C). Fixed depth also keeps GATK tractable.
  * GATK per-start downsampling is DISABLED (--max-reads-per-alignment-start 0):
    amplicon reads share a primer start, so the default cap (~50) would discard
    almost all coverage (this is why the KO pipeline also disables it).
  * All non-pindel callers are normalized with `bcftools norm -f <ref> -m -`
    (left-align + split multiallelics) before scoring, so their breakpoints line
    up with the canonical VCF-style truth anchor (simulate.make_mutant_sequence).

Outputs (under runs/sim_benchmark/exptD/<dataset>/)
  results/caller_calls.csv        one row per (condition, caller): detected,
                                  support, off-target count, runtime.
  results/caller_sensitivity.csv  per caller x grid group: TP/FP/FN and
                                  sensitivity / precision / F1 (controls supply
                                  the per-caller false-positive baseline).

Requires (in the indel_sim env): pindel/pindel2vcf, varscan, gatk (>=4),
lofreq, bcftools, samtools, pysam, pandas. Only the callers requested via
--callers need be installed.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import time

BENCH_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH_DIR)
sys.path.insert(0, os.path.join(BENCH_DIR, "utils"))

import pandas as pd  # noqa: E402

import config    # noqa: E402
import evaluate  # noqa: E402
import simulate  # noqa: E402
from align_call import bwa_runner, samtools_runner  # noqa: E402

SOURCE_EXPERIMENT = "A"
EXPERIMENT = "D"

# --- reduced comparison grid (owned here, a subset of Exp A) ------------------
CALLERS = ["pindel", "varscan", "mutect2", "lofreq"]
LENGTHS = [1, 5, 10, 20]
ALLELE_FREQS = [0.01, 0.02, 0.05]
DEPTHS = [5000, 10000]
SEEDS = [0]
# positions and read lengths are taken as-is from the source manifest.

# Raw-detectability ceiling: a call needs >= this many ALT-supporting reads.
MIN_SUPPORT = config.CEILING_E  # = 1
# Pindel operating threshold from the Exp C sweep: a fractional cutoff on
# Pindel's own AD (alt / (alt + ref)). Unlike a fixed absolute -e, a fixed 1%
# stays near-F1-optimal across depth (Exp C: 0.02@5k, 0.01@10k, 0.005@100k;
# 1% is best at 10k and close elsewhere). Applied to pindel only; the other
# callers keep their own raw output at MIN_SUPPORT.
PINDEL_OPT_FRACTION = 0.01


# --- source / destination layout --------------------------------------------

def _source_outdir(dataset):
    return os.path.join(config.OUTDIR, f"expt{SOURCE_EXPERIMENT}", dataset)


def _dest_outdir(dataset):
    return os.path.join(config.OUTDIR, f"expt{EXPERIMENT}", dataset)


def _wt_ref(source_outdir, dataset):
    return os.path.join(source_outdir, "references",
                        f"{config.amplicon_name(dataset)}_WT.fa")


def _chrom(dataset):
    return f"{config.amplicon_name(dataset)}_WT"


def _sorted_bam(cond_dir, cond_id):
    return os.path.join(cond_dir, f"{cond_id}.sorted.bam")


# --- BAM regeneration -------------------------------------------------------

def _regen_missing_bams(rows, source_outdir, dataset):
    """Re-simulate and re-align any condition whose .sorted.bam was deleted.

    Only triggered when --regen-bam is passed.  Skips conditions that already
    have a BAM.  Does NOT re-run Pindel (pindel-output VCFs still exist).
    """
    ds_cfg = config.DATASETS[dataset]
    fragment_mean = ds_cfg["fragment_mean"]
    fragment_sd   = ds_cfg["fragment_sd"]
    wt_ref = _wt_ref(source_outdir, dataset)
    refs_dir = os.path.join(source_outdir, "references")
    conds_dir = os.path.join(source_outdir, "conditions")

    missing = []
    for row in rows:
        cond_id = row["condition_id"]
        cond_dir = os.path.join(conds_dir, cond_id)
        bam = _sorted_bam(cond_dir, cond_id)
        if not os.path.exists(bam):
            missing.append(row)

    if not missing:
        print(f"[D] {dataset}: all BAMs present, nothing to regenerate.")
        return

    print(f"[D] {dataset}: regenerating {len(missing)} missing BAMs "
          f"(simulate + align, skip Pindel) ...")
    for i, row in enumerate(missing, 1):
        cond_id = row["condition_id"]
        cond_dir = os.path.join(conds_dir, cond_id)
        depth    = int(row["depth"])
        read_len = int(row["read_len"])
        seed     = int(row["seed"])
        is_ctrl  = _is_control(row)
        try:
            if is_ctrl:
                r1, r2 = simulate.simulate_control(
                    cond_dir, wt_ref, seed, depth, read_len,
                    fragment_mean, fragment_sd)
            else:
                af  = float(row["allele_freq"])
                length = int(row["length"])
                pos    = int(row["amplicon_pos"])
                indel_type = row["indel_type"]
                mut_name = (
                    f"{config.amplicon_name(dataset)}"
                    f"_{indel_type}_L{length}_p{pos}"
                )
                mut_ref = os.path.join(refs_dir, f"{mut_name}.fa")
                if not os.path.exists(mut_ref):
                    print(f"  [regen] SKIP {cond_id}: mut_ref not found: {mut_ref}")
                    continue
                r1, r2 = simulate.simulate_condition(
                    cond_dir, mut_ref, wt_ref, af, seed, depth, read_len,
                    fragment_mean, fragment_sd)
            raw_bam = bwa_runner(wt_ref, r1, r2, True)
            if not raw_bam:
                print(f"  [regen] SKIP {cond_id}: bwa failed")
                continue
            # samtools_runner short-circuits if .mpileup already exists,
            # skipping the sort/index step and leaving .sorted.bam absent.
            # Remove the stale mpileup so it re-runs the full pipeline.
            stale_mpileup = raw_bam.replace(".bam", ".mpileup")
            if os.path.exists(stale_mpileup) and not os.path.exists(
                    raw_bam.replace(".bam", ".sorted.bam")):
                os.remove(stale_mpileup)
            _, sorted_bam, _ = samtools_runner(raw_bam, wt_ref)
            if not sorted_bam:
                print(f"  [regen] SKIP {cond_id}: samtools failed")
                continue
            # Clean up FASTQs to save disk.
            for fq in [r1, r2]:
                try:
                    os.remove(fq)
                except OSError:
                    pass
            print(f"  [regen] {i}/{len(missing)} done: {cond_id}")
        except Exception as exc:  # noqa: BLE001
            print(f"  [regen] ERROR {cond_id}: {exc}")
    print(f"[D] {dataset}: BAM regeneration complete.")


# --- source manifest selection ----------------------------------------------

def _is_control(row):
    return str(row.get("is_control", "")).strip() in ("1", "True", "true")


def _select_rows(manifest):
    """Keep grid-matching mutants and the controls that back their precision."""
    mutants, controls = [], []
    for row in manifest:
        read_len = int(row["read_len"])
        depth = int(row["depth"])
        seed = int(row["seed"])
        if depth not in DEPTHS or seed not in SEEDS:
            continue
        if _is_control(row):
            controls.append(row)
        elif (int(row["length"]) in LENGTHS
              and float(row["allele_freq"]) in ALLELE_FREQS):
            mutants.append(row)
    return mutants, controls


# --- generic normalized-VCF parsing -----------------------------------------

def _bcftools_norm(raw_vcf, ref, out_vcf):
    """Left-align + split multiallelics so breakpoints match the truth anchor.

    Some callers (e.g. VarScan) omit ``##contig`` header lines, which makes
    ``bcftools norm`` fail. Reheader from the reference ``.fai`` first to inject
    them; this is a no-op when the VCF already declares the contig.
    """
    reheader = subprocess.Popen(
        ["bcftools", "reheader", "--fai", ref + ".fai", raw_vcf],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    with open(out_vcf, "w") as out:
        subprocess.check_call(
            ["bcftools", "norm", "-f", ref, "-m", "-"],
            stdin=reheader.stdout, stdout=out, stderr=subprocess.DEVNULL)
    reheader.stdout.close()
    if reheader.wait() != 0:
        raise subprocess.CalledProcessError(reheader.returncode, "bcftools reheader")


def _parse_norm_vcf(path, support_kind):
    """Parse a normalized VCF into indel call dicts {pos, svtype, abs_svlen, support}."""
    calls = []
    if not os.path.exists(path):
        return calls
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            ref, alt = f[3], f[4]
            if alt in (".", "") or "," in alt:
                continue
            dlen = len(alt) - len(ref)
            if dlen == 0:
                continue
            calls.append({
                "pos": int(f[1]),
                "svtype": "INS" if dlen > 0 else "DEL",
                "abs_svlen": abs(dlen),
                "support": _extract_support(f, support_kind),
            })
    return calls


def _extract_support(f, kind):
    """ALT-supporting read count, per-caller field layout."""
    try:
        if kind == "varscan_ad":            # FORMAT AD = ALT reads
            d = dict(zip(f[8].split(":"), f[9].split(":")))
            return int(d.get("AD", "0"))
        if kind == "gatk_ad":               # FORMAT AD = ref,alt
            d = dict(zip(f[8].split(":"), f[9].split(":")))
            return int(d.get("AD", "0,0").split(",")[1])
        if kind == "lofreq":                # INFO DP4 = rf,rr,af,ar (fallback AF*DP)
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            if "DP4" in info:
                p = info["DP4"].split(",")
                return int(p[2]) + int(p[3])
            if "AF" in info and "DP" in info:
                return int(round(float(info["AF"]) * int(info["DP"])))
            return 0
    except (ValueError, IndexError):
        return 0
    return 0


# --- read-group helper (Mutect2 requires a sample) --------------------------

def _ensure_rg_bam(bam, workdir):
    """Return a BAM guaranteed to carry an @RG (adds one only if missing)."""
    header = subprocess.check_output(
        ["samtools", "view", "-H", bam], text=True)
    if any(line.startswith("@RG") for line in header.splitlines()):
        return bam
    out = os.path.join(workdir, "rg.bam")
    subprocess.check_call(
        ["samtools", "addreplacerg",
         "-r", "@RG\tID:sim\tSM:sim\tPL:ILLUMINA",
         "-o", out, bam])
    subprocess.check_call(["samtools", "index", out])
    return out


# --- per-caller runners ------------------------------------------------------

def _run_pindel(cond_dir, bam, ref, chrom, workdir):
    """Re-use Exp A's base -e 1 pindel-output VCFs (scoring applies a fractional
    AD cutoff, PINDEL_OPT_FRACTION). Carries ref_reads so alt / (alt + ref) can
    be computed at score time."""
    pindel_dir = os.path.join(cond_dir, "pindel-output")
    parsed = evaluate.parse_condition(pindel_dir)
    calls = []
    for svtype in ("DEL", "INS"):
        for c in parsed[svtype]:
            calls.append({
                "pos": c["pos"], "svtype": c["svtype"],
                "abs_svlen": c["abs_svlen"], "support": c["alt_reads"],
                "ref_reads": c["ref_reads"],
            })
    return calls


def _run_varscan(cond_dir, bam, ref, chrom, workdir):
    mpileup = os.path.join(workdir, "cond.mpileup")
    with open(mpileup, "w") as out:
        subprocess.check_call(
            ["samtools", "mpileup", "-B", "-A", "-d", "1000000",
             "-f", ref, "-r", chrom, bam],
            stdout=out, stderr=subprocess.DEVNULL)
    raw = os.path.join(workdir, "varscan_indel.vcf")
    with open(raw, "w") as out:
        subprocess.check_call(
            ["varscan", "mpileup2indel", mpileup,
             "--min-var-freq", "0.001", "--min-reads2", "1",
             "--min-avg-qual", "0", "--p-value", "0.99",
             "--strand-filter", "0", "--output-vcf", "1"],
            stdout=out, stderr=subprocess.DEVNULL)
    norm = os.path.join(workdir, "varscan.norm.vcf")
    _bcftools_norm(raw, ref, norm)
    return _parse_norm_vcf(norm, "varscan_ad")


def _run_mutect2(cond_dir, bam, ref, chrom, workdir):
    ref_dict = os.path.splitext(ref)[0] + ".dict"
    if not os.path.exists(ref_dict) and not os.path.exists(ref + ".dict"):
        subprocess.check_call(
            ["gatk", "CreateSequenceDictionary", "-R", ref],
            stderr=subprocess.DEVNULL)
    rg_bam = _ensure_rg_bam(bam, workdir)
    raw = os.path.join(workdir, "mutect2.vcf.gz")
    subprocess.check_call(
        ["gatk", "Mutect2", "-R", ref, "-I", rg_bam, "-L", chrom,
         "--max-reads-per-alignment-start", "0", "-O", raw],
        stderr=subprocess.DEVNULL)
    norm = os.path.join(workdir, "mutect2.norm.vcf")
    _bcftools_norm(raw, ref, norm)
    return _parse_norm_vcf(norm, "gatk_ad")


def _run_lofreq(cond_dir, bam, ref, chrom, workdir):
    rg_bam = _ensure_rg_bam(bam, workdir)
    iq_bam = os.path.join(workdir, "iq.bam")
    subprocess.check_call(
        ["lofreq", "indelqual", "--dindel", "-f", ref, "-o", iq_bam, rg_bam],
        stderr=subprocess.DEVNULL)
    subprocess.check_call(["samtools", "index", iq_bam])
    raw = os.path.join(workdir, "lofreq.vcf")
    subprocess.check_call(
        ["lofreq", "call", "--call-indels", "-f", ref, "-r", chrom,
         "-o", raw, iq_bam],
        stderr=subprocess.DEVNULL)
    norm = os.path.join(workdir, "lofreq.norm.vcf")
    _bcftools_norm(raw, ref, norm)
    return _parse_norm_vcf(norm, "lofreq")


_RUNNERS = {
    "pindel": _run_pindel,
    "varscan": _run_varscan,
    "mutect2": _run_mutect2,
    "lofreq": _run_lofreq,
}


# --- caching -----------------------------------------------------------------

def _cache_paths(dest_outdir, cond_id, caller):
    d = os.path.join(dest_outdir, "conditions", cond_id)
    return (os.path.join(d, f"{caller}.calls.tsv"),
            os.path.join(d, f"{caller}.runtime"), d)


def _load_calls(calls_tsv, runtime_file):
    # ref_reads (5th column) is only written for pindel; absent for callers that
    # are scored on absolute support, where it is left as "".
    calls = []
    with open(calls_tsv) as fh:
        next(fh, None)  # header
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            ref_reads = int(p[4]) if len(p) >= 5 and p[4] != "" else ""
            calls.append({"pos": int(p[0]), "svtype": p[1],
                          "abs_svlen": int(p[2]), "support": int(p[3]),
                          "ref_reads": ref_reads})
    runtime = ""
    if os.path.exists(runtime_file):
        runtime = open(runtime_file).read().strip()
    return calls, runtime


def _save_calls(calls_tsv, runtime_file, calls, runtime):
    with open(calls_tsv, "w") as fh:
        fh.write("pos\tsvtype\tabs_svlen\tsupport\tref_reads\n")
        for c in calls:
            fh.write(f"{c['pos']}\t{c['svtype']}\t{c['abs_svlen']}\t"
                     f"{c['support']}\t{c.get('ref_reads', '')}\n")
    with open(runtime_file, "w") as fh:
        fh.write(f"{runtime:.3f}")


def _caller_calls(caller, cond_dir, bam, ref, chrom, dest_outdir, cond_id,
                  force=False):
    """Return (calls, runtime_s), running the caller or loading the cache."""
    calls_tsv, runtime_file, work_parent = _cache_paths(dest_outdir, cond_id, caller)
    if os.path.exists(calls_tsv) and not force:
        return _load_calls(calls_tsv, runtime_file)
    workdir = os.path.join(work_parent, caller)
    os.makedirs(workdir, exist_ok=True)
    start = time.time()
    calls = _RUNNERS[caller](cond_dir, bam, ref, chrom, workdir)
    runtime = time.time() - start
    _save_calls(calls_tsv, runtime_file, calls, runtime)
    return calls, f"{runtime:.3f}"


# --- scoring -----------------------------------------------------------------

def _passes(call, min_support, min_fraction=None):
    """Keep a call if it clears the operating threshold.

    When min_fraction is set (pindel), filter on the AD fraction
    alt / (alt + ref); otherwise filter on absolute alt support.
    """
    if min_fraction is not None:
        ref = call.get("ref_reads", "")
        ref = int(ref) if ref != "" else 0
        denom = call["support"] + ref
        return denom > 0 and (call["support"] / denom) >= min_fraction
    return call["support"] >= min_support


def _score_mutant(calls, expected_svtype, expected_len, truth_pos, min_support,
                  min_fraction=None):
    kept = [c for c in calls if _passes(c, min_support, min_fraction)]
    exact = [c for c in kept
             if evaluate.is_exact(c, expected_svtype, expected_len, truth_pos)]
    same_type = [c for c in kept if c["svtype"] == expected_svtype]
    offtarget = [c for c in same_type if c not in exact]
    support = max((c["support"] for c in exact), default="")
    return int(bool(exact)), len(offtarget), support


def _score_control(calls, min_support, min_fraction=None):
    kept = [c for c in calls if _passes(c, min_support, min_fraction)]
    return {"DEL": sum(1 for c in kept if c["svtype"] == "DEL"),
            "INS": sum(1 for c in kept if c["svtype"] == "INS")}


# --- aggregation -------------------------------------------------------------

def _confusion(call_rows, control_fp):
    """Per (caller x grid group) TP/FP/FN + sensitivity/precision/F1."""
    df = pd.DataFrame([r for r in call_rows if r["is_control"] == 0])
    if df.empty:
        return pd.DataFrame()
    group_cols = ["caller", "dataset", "indel_type", "length", "allele_freq",
                  "read_len", "depth"]
    out = []
    for keys, g in df.groupby(group_cols):
        caller, dataset, indel_type, length, af, read_len, depth = keys
        svtype = "DEL" if indel_type == "del" else "INS"
        tp = int(g["detected"].sum())
        n = len(g)
        fn = n - tp
        fp = int(g["n_offtarget"].sum())
        fp += control_fp.get((caller, dataset, read_len, depth), {}).get(svtype, 0)
        sens = evaluate._safe_div(tp, tp + fn)
        prec = evaluate._safe_div(tp, tp + fp)
        f1 = ""
        if prec != "" and sens != "" and (prec + sens) > 0:
            f1 = round(2 * prec * sens / (prec + sens), 4)
        out.append({
            "caller": caller, "dataset": dataset, "indel_type": indel_type,
            "length": length, "allele_freq": af, "read_len": read_len,
            "depth": depth, "n": n, "TP": tp, "FP": fp, "FN": fn,
            "sensitivity": round(sens, 4) if sens != "" else "",
            "precision": round(prec, 4) if prec != "" else "",
            "f1": f1,
        })
    return pd.DataFrame(out).sort_values(group_cols)


# --- per-dataset driver ------------------------------------------------------

def _score_dataset(dataset, callers, force=False, regen_bam=False):
    source_outdir = _source_outdir(dataset)
    dest_outdir = _dest_outdir(dataset)
    manifest = evaluate.load_manifest(source_outdir, SOURCE_EXPERIMENT)
    if manifest is None:
        print(f"[D] {dataset}: no Exp {SOURCE_EXPERIMENT} manifest; skipping.")
        return None

    mutants, controls = _select_rows(manifest)

    if regen_bam:
        _regen_missing_bams(mutants + controls, source_outdir, dataset)
    ref = _wt_ref(source_outdir, dataset)
    chrom = _chrom(dataset)
    os.makedirs(os.path.join(dest_outdir, "results"), exist_ok=True)
    print(f"[D] {dataset}: {len(mutants)} mutant + {len(controls)} control "
          f"conditions x {len(callers)} callers")

    call_rows = []
    control_fp = {}
    total = len(mutants) + len(controls)
    done = 0

    for row in mutants:
        cond_id = row["condition_id"]
        cond_dir = os.path.join(source_outdir, "conditions", cond_id)
        bam = _sorted_bam(cond_dir, cond_id)
        expected_svtype = row["expected_svtype"]
        length = int(row["length"])
        truth_pos = int(row["anchor_pos"])
        for caller in callers:
            calls, runtime = _caller_calls(caller, cond_dir, bam, ref, chrom,
                                           dest_outdir, cond_id, force=force)
            min_fraction = PINDEL_OPT_FRACTION if caller == "pindel" else None
            detected, n_off, support = _score_mutant(
                calls, expected_svtype, length, truth_pos, MIN_SUPPORT,
                min_fraction)
            call_rows.append({
                "condition_id": cond_id, "caller": caller, "dataset": dataset,
                "is_control": 0, "indel_type": row["indel_type"], "length": length,
                "allele_freq": float(row["allele_freq"]),
                "amplicon_pos": int(row["amplicon_pos"]),
                "read_len": int(row["read_len"]), "depth": int(row["depth"]),
                "seed": int(row["seed"]), "detected": detected,
                "n_offtarget": n_off, "support": support, "runtime_s": runtime,
            })
        done += 1
        if done % 25 == 0:
            print(f"[D]   {done}/{total} conditions scored")

    for row in controls:
        cond_id = row["condition_id"]
        cond_dir = os.path.join(source_outdir, "conditions", cond_id)
        bam = _sorted_bam(cond_dir, cond_id)
        read_len, depth = int(row["read_len"]), int(row["depth"])
        for caller in callers:
            calls, runtime = _caller_calls(caller, cond_dir, bam, ref, chrom,
                                           dest_outdir, cond_id, force=force)
            min_fraction = PINDEL_OPT_FRACTION if caller == "pindel" else None
            fp = _score_control(calls, MIN_SUPPORT, min_fraction)
            control_fp.setdefault((caller, dataset, read_len, depth), {}).update(fp)
            call_rows.append({
                "condition_id": cond_id, "caller": caller, "dataset": dataset,
                "is_control": 1, "indel_type": "", "length": "",
                "allele_freq": 0.0, "amplicon_pos": "",
                "read_len": read_len, "depth": depth, "seed": int(row["seed"]),
                "detected": "", "n_offtarget": fp["DEL"] + fp["INS"],
                "support": "", "runtime_s": runtime,
            })
        done += 1

    calls_df = pd.DataFrame(call_rows).sort_values(["condition_id", "caller"])
    calls_path = os.path.join(dest_outdir, "results", "caller_calls.csv")
    calls_df.to_csv(calls_path, index=False)

    conf = _confusion(call_rows, control_fp)
    conf_path = os.path.join(dest_outdir, "results", "caller_sensitivity.csv")
    conf.to_csv(conf_path, index=False)
    print(f"[D] {dataset}: wrote {calls_path}\n[D] {dataset}: wrote {conf_path}")
    return conf


def main(datasets=None, callers=None, force=False, regen_bam=False):
    if datasets is None:
        datasets = list(config.DATASETS)
    if callers is None:
        callers = CALLERS
    for dataset in datasets:
        _score_dataset(dataset, callers, force=force, regen_bam=regen_bam)
    print("[D] all datasets done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Experiment D: indel-caller sensitivity comparison "
                    "(re-scores Exp A BAMs; no re-simulation).")
    parser.add_argument("-d", "--dataset", choices=list(config.DATASETS),
                        help="Run one dataset only (omit to run all).")
    parser.add_argument("--callers", default=",".join(CALLERS),
                        help=f"Comma-separated subset of {CALLERS}.")
    parser.add_argument("--force", action="store_true",
                        help="Re-run callers even if a cached calls.tsv exists.")
    parser.add_argument("--regen-bam", action="store_true",
                        help="Re-simulate + re-align any condition whose "
                             ".sorted.bam was deleted (e.g. after disk cleanup). "
                             "Skips Pindel (pindel-output VCFs must already exist).")
    args = parser.parse_args()
    chosen = [c.strip() for c in args.callers.split(",") if c.strip()]
    unknown = [c for c in chosen if c not in _RUNNERS]
    if unknown:
        parser.error(f"unknown caller(s): {unknown}; choose from {CALLERS}")
    main(
        datasets=[args.dataset] if args.dataset else None,
        callers=chosen,
        force=args.force,
        regen_bam=args.regen_bam,
    )
