#!/usr/bin/env python3
"""Summarize the output of compare_runs.py (or compare_depth_sweep.py).

Reads a comparison output directory produced by compare_runs.py and prints a
human-readable summary: per-sample matched/old-only/new-only counts for SNPs
and indels, plus any notable VarFreq deltas in matched calls.

Usage:
  python comparison_summary.py comparison/AvB1_HC/

Reconstructed Aug 2026 from shell history after accidental deletion.
"""
import csv
import glob
import os
import sys


NOTABLE_DELTA_PP = 2.0  # flag matched calls with |delta| >= 2 percentage points


def _read_tsv(path):
    """Return list of dicts from a TSV file."""
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _parse_delta(s):
    """'+1.23pp' -> 1.23"""
    try:
        return float(s.replace("pp", "").strip())
    except (ValueError, TypeError, AttributeError):
        return None


def summarize(comparison_dir):
    """Print a summary of all comparison TSVs in *comparison_dir*."""
    summary_path = os.path.join(comparison_dir, "summary.tsv")
    if os.path.isfile(summary_path):
        rows = _read_tsv(summary_path)
        if rows:
            print(f"{'Sample':<40s} {'SNP':>8s} {'IND':>8s}  "
                  f"{'S.match':>7s} {'S.old':>6s} {'S.new':>6s}  "
                  f"{'I.match':>7s} {'I.old':>6s} {'I.new':>6s}")
            print("-" * 100)
            for r in rows:
                print(f"{r.get('sample','?'):<40s} "
                      f"{r.get('snp_old',''):>8s} {r.get('indel_old',''):>8s}  "
                      f"{r.get('snp_matched',''):>7s} "
                      f"{r.get('snp_old_only',''):>6s} "
                      f"{r.get('snp_new_only',''):>6s}  "
                      f"{r.get('indel_matched',''):>7s} "
                      f"{r.get('indel_old_only',''):>6s} "
                      f"{r.get('indel_new_only',''):>6s}")
            print()

    # Look for notable frequency deltas in *.matched.tsv files.
    matched_files = sorted(glob.glob(os.path.join(comparison_dir, "*.matched.tsv")))
    notable = []
    for mf in matched_files:
        basename = os.path.basename(mf)
        for row in _read_tsv(mf):
            delta = _parse_delta(row.get("DeltaFreq", ""))
            if delta is not None and abs(delta) >= NOTABLE_DELTA_PP:
                notable.append((basename, row))

    if notable:
        print(f"Notable VarFreq shifts (|delta| >= {NOTABLE_DELTA_PP}pp):")
        for fname, row in notable:
            print(f"  {fname}  pos {row.get('Position','')}  "
                  f"{row.get('Ref','')}->{row.get('VarAllele','')}  "
                  f"old={row.get('OldFreq','')}  new={row.get('NewFreq','')}  "
                  f"delta={row.get('DeltaFreq','')}")
        print()
    elif matched_files:
        print(f"No matched calls with |VarFreq delta| >= {NOTABLE_DELTA_PP}pp.\n")

    # Report any old-only or new-only files present.
    old_only_files = sorted(glob.glob(os.path.join(comparison_dir, "*.old_only.tsv")))
    new_only_files = sorted(glob.glob(os.path.join(comparison_dir, "*.new_only.tsv")))
    if old_only_files:
        print("Old-only calls (variants lost in new run):")
        for f in old_only_files:
            rows = _read_tsv(f)
            print(f"  {os.path.basename(f)}: {len(rows)} variant(s)")
            for r in rows[:5]:
                print(f"    pos {r.get('Position','')}  "
                      f"{r.get('Ref','')}->{r.get('VarAllele','')}  "
                      f"freq={r.get('VarFreq','')}")
            if len(rows) > 5:
                print(f"    ... and {len(rows)-5} more")
        print()

    if new_only_files:
        print("New-only calls (variants gained in new run):")
        for f in new_only_files:
            rows = _read_tsv(f)
            print(f"  {os.path.basename(f)}: {len(rows)} variant(s)")
            for r in rows[:5]:
                print(f"    pos {r.get('Position','')}  "
                      f"{r.get('Ref','')}->{r.get('VarAllele','')}  "
                      f"freq={r.get('VarFreq','')}")
            if len(rows) > 5:
                print(f"    ... and {len(rows)-5} more")
        print()


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <comparison_output_dir>", file=sys.stderr)
        return 1
    comparison_dir = sys.argv[1]
    if not os.path.isdir(comparison_dir):
        print(f"Not a directory: {comparison_dir}", file=sys.stderr)
        return 1
    summarize(comparison_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
