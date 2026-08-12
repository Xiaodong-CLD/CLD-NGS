#!/usr/bin/env python3
# 2026-06-varscan-interface/compare_varscan_interface.py
# Diagnostic: does switching VarScan's legacy `pileup2*` callers to the modern
# `mpileup2*` callers change SNV / indel calls?
#
# For each sample's existing d<D>.mpileup (from the 2026-06-mpileup-depth-cap
# sweep), runs both interfaces with IDENTICAL explicit parameters (so only the
# interface differs, not its different built-in defaults), normalizes the two different
# output formats to a common {position: (ref, var, varfreq)} call set, and
# reports positions only in one interface plus VarFreq differences on shared
# positions. SNPs and indels are compared separately.
#
# Important context (why identical params are forced):
#   - pileup2snp default --min-var-freq is 0.01 (1%); mpileup2snp default is 0.2
#     (20%) -> without matching, mpileup2snp drops all low-freq amplicon calls.
#   - The two interfaces emit DIFFERENT output formats:
#       pileup2*  -> CONS format: Chrom Pos Ref Cons Reads1 Reads2 VarFreq ... VarAllele(col18)
#       mpileup2* -> packed format: Chrom Pos Ref Var Cons:Cov:Reads1:Reads2:Freq:P-value ...
#
# Usage:
#   conda activate variant-calling
#   python compare_varscan_interface.py <repo>/runs/T3_rebuild <repo>/runs/T4_varscan_update
#   python compare_varscan_interface.py <repo>/runs/T3_rebuild <repo>/runs/T4_varscan_update --depths 0 8000

import argparse
import glob
import os
import subprocess
import sys

# Parameters forced identical on BOTH interfaces, isolating the interface from
# its differing defaults. SNP keeps the production --min-coverage 1500; indel
# uses VarScan's default min-coverage 8 (matches the production pileup2indel call).
COMMON_SNP = ["--min-coverage", "1500", "--min-var-freq", "0.01",
              "--min-reads2", "2", "--strand-filter", "0"]
COMMON_INDEL = ["--min-coverage", "8", "--min-var-freq", "0.01",
                "--min-reads2", "2", "--strand-filter", "0"]
VARFREQ_REPORT_THRESHOLD = 2.0  # percent; also summarize calls at the pipeline's 2% cutoff


def run_varscan(subcommand, mpileup_path, params, out_path):
    """Run `varscan <subcommand> <mpileup> <params>` into out_path (cached)."""
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return out_path
    with open(out_path, "w") as out:
        subprocess.check_call(["varscan", subcommand, mpileup_path] + params, stdout=out)
    return out_path


def parse_freq(token):
    """'2.34%' -> 2.34 (float). Returns None if unparseable."""
    token = token.strip().rstrip("%")
    try:
        return float(token)
    except ValueError:
        return None


def parse_pileup2_format(path):
    """Parse pileup2snp/pileup2indel (CONS) output.

    Columns: 0 Chrom, 1 Position, 2 Ref, 3 Cons, 4 Reads1, 5 Reads2,
             6 VarFreq, ..., 18 VarAllele.
    Returns {position(int): (ref, var_allele, varfreq_float)}.
    """
    calls = {}
    with open(path) as fh:
        next(fh, None)  # header
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 19 or not cols[1].isdigit():
                continue
            freq = parse_freq(cols[6])
            if freq is None:
                continue
            calls[int(cols[1])] = (cols[2], cols[18], freq)
    return calls


def parse_mpileup2_format(path):
    """Parse mpileup2snp/mpileup2indel (packed) output.

    Columns: 0 Chrom, 1 Position, 2 Ref, 3 Var,
             4 Cons:Cov:Reads1:Reads2:Freq:P-value (aggregate), ...
    Returns {position(int): (ref, var, varfreq_float)}.
    """
    calls = {}
    with open(path) as fh:
        next(fh, None)  # header
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5 or not cols[1].isdigit():
                continue
            agg = cols[4].split(":")
            if len(agg) < 5:
                continue
            freq = parse_freq(agg[4])
            if freq is None:
                continue
            calls[int(cols[1])] = (cols[2], cols[3], freq)
    return calls


def compare_calls(pileup_calls, mpileup_calls):
    """Compare two {pos: (ref, var, freq)} dicts. Returns a summary dict."""
    pset, mset = set(pileup_calls), set(mpileup_calls)
    shared = pset & mset
    max_freq_diff = 0.0
    var_mismatch = 0
    for pos in shared:
        _, pv, pf = pileup_calls[pos]
        _, mv, mf = mpileup_calls[pos]
        max_freq_diff = max(max_freq_diff, abs(pf - mf))
        if pv != mv:
            var_mismatch += 1
    p_thr = {pos for pos in pset if pileup_calls[pos][2] >= VARFREQ_REPORT_THRESHOLD}
    m_thr = {pos for pos in mset if mpileup_calls[pos][2] >= VARFREQ_REPORT_THRESHOLD}
    return {
        "n_pileup": len(pset),
        "n_mpileup": len(mset),
        "n_shared": len(shared),
        "only_pileup": sorted(pset - mset),
        "only_mpileup": sorted(mset - pset),
        "var_mismatch": var_mismatch,
        "max_freq_diff": max_freq_diff,
        "n_pileup_ge_thr": len(p_thr),
        "n_mpileup_ge_thr": len(m_thr),
        "thr_only_pileup": sorted(p_thr - m_thr),
        "thr_only_mpileup": sorted(m_thr - p_thr),
        "identical": pset == mset and var_mismatch == 0 and max_freq_diff == 0.0,
    }


def discover_mpileups(run_dir, depths):
    """Yield (sample, depth, mpileup_path) for existing d<D>.mpileup files."""
    for sweep_dir in sorted(glob.glob(os.path.join(run_dir, "Results", "*", "depth-sweep"))):
        sample = os.path.basename(os.path.dirname(sweep_dir))
        for d in depths:
            mp = os.path.join(sweep_dir, "{}.d{}.mpileup".format(sample, d))
            if os.path.isfile(mp):
                yield sample, d, mp


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", help="Run directory containing Results/*/depth-sweep/")
    ap.add_argument("out_dir", help="Directory to write comparison outputs into")
    ap.add_argument("--depths", nargs="+", type=int, default=[0, 8000],
                    help="Depth caps to compare (default: 0 8000)")
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    report_path = os.path.join(args.out_dir, "interface_comparison_report.txt")
    overall_identical = True
    lines = []

    def emit(text=""):
        lines.append(text)
        print(text)

    emit("VarScan interface comparison: pileup2* (legacy) vs mpileup2* (modern)")
    emit("Run: {}   Depths: {}".format(args.run_dir, args.depths))
    emit("Forced-identical SNP params:   {}".format(" ".join(COMMON_SNP)))
    emit("Forced-identical indel params: {}".format(" ".join(COMMON_INDEL)))
    emit("=" * 78)

    work = list(discover_mpileups(args.run_dir, args.depths))
    if not work:
        sys.stderr.write("No d<D>.mpileup files found under {}/Results/*/depth-sweep/\n".format(
            args.run_dir))
        return 1

    for sample, d, mp in work:
        sample_out = os.path.join(args.out_dir, sample)
        os.makedirs(sample_out, exist_ok=True)
        emit("")
        emit("Sample {}  |  d={}".format(sample, d))
        emit("-" * 78)

        for vartype, legacy, modern, params, plegacy, pmodern in (
            ("SNP", "pileup2snp", "mpileup2snp", COMMON_SNP,
             parse_pileup2_format, parse_mpileup2_format),
            ("INDEL", "pileup2indel", "mpileup2indel", COMMON_INDEL,
             parse_pileup2_format, parse_mpileup2_format),
        ):
            legacy_out = os.path.join(sample_out, "{}.d{}.{}.txt".format(sample, d, legacy))
            modern_out = os.path.join(sample_out, "{}.d{}.{}.txt".format(sample, d, modern))
            run_varscan(legacy, mp, params, legacy_out)
            run_varscan(modern, mp, params, modern_out)
            cmp = compare_calls(plegacy(legacy_out), pmodern(modern_out))
            overall_identical = overall_identical and cmp["identical"]

            emit("  {}: {} reported {} calls, {} reported {} calls; shared {}".format(
                vartype, legacy, cmp["n_pileup"], modern, cmp["n_mpileup"], cmp["n_shared"]))
            emit("    only in {}: {}".format(legacy, cmp["only_pileup"] or "none"))
            emit("    only in {}: {}".format(modern, cmp["only_mpileup"] or "none"))
            emit("    var-allele mismatches on shared: {};  max |VarFreq diff|: {:.4f} pp".format(
                cmp["var_mismatch"], cmp["max_freq_diff"]))
            emit("    at >= {:.1f}%: {} {} vs {} {}  (only-{}: {}, only-{}: {})".format(
                VARFREQ_REPORT_THRESHOLD, cmp["n_pileup_ge_thr"], legacy,
                cmp["n_mpileup_ge_thr"], modern,
                legacy, cmp["thr_only_pileup"] or "none",
                modern, cmp["thr_only_mpileup"] or "none"))
            emit("    => {}".format("IDENTICAL calls" if cmp["identical"]
                                    else "DIFFERENT calls"))

    emit("")
    emit("=" * 78)
    emit("OVERALL: {}".format(
        "pileup2* and mpileup2* produce IDENTICAL calls under matched params."
        if overall_identical else
        "pileup2* and mpileup2* DIFFER. Do NOT switch the wrapper without "
        "handling format/parameter differences."))

    with open(report_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n[written] {}".format(report_path))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
