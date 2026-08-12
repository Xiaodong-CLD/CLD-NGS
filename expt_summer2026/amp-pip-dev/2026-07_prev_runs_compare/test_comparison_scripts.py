#!/usr/bin/env python3
"""Smoke test for the reconstructed 2026-07_prev_runs_compare scripts.

Creates synthetic amp-pip-dev Results/ layouts with VarScan SNP/indel/consensus
files, runs compare_runs.py and comparison_summary.py against them, and verifies
outputs are correct. Also tests compare_depth_sweep.py with synthetic depth-sweep
consensus files.

Run from the 2026-07_prev_runs_compare/ directory:
    python test_comparison_scripts.py

Or with pytest:
    pytest test_comparison_scripts.py -v
"""
import csv
import os
import shutil
import subprocess
import sys
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


# ── Helpers ───────────────────────────────────────────────────────────────────

def _write_varscan_snp(path, rows):
    """Write a VarScan pileup2snp-style TSV."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    header = ["Chrom", "Position", "Ref", "Cons", "Reads1", "Reads2",
              "VarFreq", "Strands1", "Strands2", "Qual1", "Qual2",
              "Pvalue", "MapQual1", "MapQual2", "Reads1Plus", "Reads1Minus",
              "Reads2Plus", "Reads2Minus", "VarAllele"]
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            line = [
                r.get("Chrom", "ref"), str(r["Position"]), r["Ref"],
                r.get("Cons", r["VarAllele"]),
                str(r.get("Reads1", 10000)), str(r.get("Reads2", 500)),
                r["VarFreq"],
                "2", "2", "40", "40", "0.001", "40", "40",
                "5000", "5000", "250", "250", r["VarAllele"],
            ]
            fh.write("\t".join(line) + "\n")


def _write_varscan_indel(path, rows):
    """Write a VarScan pileup2indel-style TSV."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    header = ["Chrom", "Position", "Ref", "Cons", "Reads1", "Reads2",
              "VarFreq", "Strands1", "Strands2", "Qual1", "Qual2",
              "Pvalue", "MapQual1", "MapQual2", "Reads1Plus", "Reads1Minus",
              "Reads2Plus", "Reads2Minus", "VarAllele"]
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            line = [
                r.get("Chrom", "ref"), str(r["Position"]), r["Ref"],
                r.get("Cons", r["VarAllele"]),
                str(r.get("Reads1", 10000)), str(r.get("Reads2", 800)),
                r["VarFreq"],
                "2", "2", "40", "40", "0.001", "40", "40",
                "5000", "5000", "400", "400", r["VarAllele"],
            ]
            fh.write("\t".join(line) + "\n")


def _write_varscan_consensus(path, rows):
    """Write a VarScan pileup2cns-style TSV."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    header = ["Chrom", "Position", "Ref", "Cons", "Reads1", "Reads2", "VarFreq",
              "Strands1", "Strands2", "Qual1", "Qual2", "Pvalue"]
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            line = [
                r.get("Chrom", "ref"), str(r["Position"]), r["Ref"],
                r.get("Cons", r.get("Var", "A")),
                str(r.get("Reads1", 10000)), str(r.get("Reads2", 200)),
                r["VarFreq"],
                "2", "2", "40", "40", "0.001",
            ]
            fh.write("\t".join(line) + "\n")


def _build_results_dir(base, sample, snp_rows, indel_rows):
    """Build a Results/<sample>/varscan-output/ tree with SNP + indel files."""
    vdir = os.path.join(base, sample, "varscan-output")
    os.makedirs(vdir, exist_ok=True)
    _write_varscan_snp(os.path.join(vdir, f"{sample}.snp.vcf"), snp_rows)
    _write_varscan_indel(os.path.join(vdir, f"{sample}.indel.vcf"), indel_rows)


def _read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ── Test: compare_runs.py ─────────────────────────────────────────────────────

def test_compare_runs():
    """End-to-end test: create two Results dirs, run compare_runs.py, check output."""
    tmp = tempfile.mkdtemp(prefix="test_compare_runs_")
    try:
        old_dir = os.path.join(tmp, "old_results")
        new_dir = os.path.join(tmp, "new_results")
        out_dir = os.path.join(tmp, "comparison")

        # Sample: SampleA — shared SNPs (one with freq delta), one old-only, one new-only
        _build_results_dir(old_dir, "SampleA",
            snp_rows=[
                {"Position": 100, "Ref": "A", "VarAllele": "G", "VarFreq": "5.0%"},
                {"Position": 200, "Ref": "C", "VarAllele": "T", "VarFreq": "10.0%"},
                {"Position": 300, "Ref": "G", "VarAllele": "A", "VarFreq": "2.5%"},
            ],
            indel_rows=[
                {"Position": 50, "Ref": "A", "VarAllele": "-AAA", "VarFreq": "8.0%"},
            ])

        _build_results_dir(new_dir, "SampleA",
            snp_rows=[
                {"Position": 100, "Ref": "A", "VarAllele": "G", "VarFreq": "5.2%"},
                {"Position": 200, "Ref": "C", "VarAllele": "T", "VarFreq": "13.0%"},
                {"Position": 400, "Ref": "T", "VarAllele": "C", "VarFreq": "1.8%"},
            ],
            indel_rows=[
                {"Position": 50, "Ref": "A", "VarAllele": "-AAA", "VarFreq": "7.5%"},
            ])

        # Run compare_runs.py
        result = subprocess.run(
            [sys.executable, os.path.join(SCRIPT_DIR, "compare_runs.py"),
             "--old-results", old_dir, "--new-results", new_dir,
             "--output-dir", out_dir],
            capture_output=True, text=True)
        assert result.returncode == 0, f"compare_runs.py failed:\n{result.stderr}"

        # Verify summary.tsv
        summary = _read_tsv(os.path.join(out_dir, "summary.tsv"))
        assert len(summary) == 1
        row = summary[0]
        assert row["sample"] == "SampleA"
        assert row["snp_matched"] == "2"   # pos 100, 200
        assert row["snp_old_only"] == "1"  # pos 300
        assert row["snp_new_only"] == "1"  # pos 400
        assert row["indel_matched"] == "1"  # pos 50
        assert row["indel_old_only"] == "0"
        assert row["indel_new_only"] == "0"

        # Verify matched file exists and has freq delta
        matched = _read_tsv(os.path.join(out_dir, "SampleA.snp.matched.tsv"))
        assert len(matched) == 2
        pos200 = [r for r in matched if r["Position"] == "200"][0]
        assert "3.00" in pos200["DeltaFreq"]  # 13.0 - 10.0 = +3.00pp

        # Verify old_only / new_only files
        old_only = _read_tsv(os.path.join(out_dir, "SampleA.snp.old_only.tsv"))
        assert len(old_only) == 1
        assert old_only[0]["Position"] == "300"

        new_only = _read_tsv(os.path.join(out_dir, "SampleA.snp.new_only.tsv"))
        assert len(new_only) == 1
        assert new_only[0]["Position"] == "400"

        print("PASS: test_compare_runs")

    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ── Test: comparison_summary.py ───────────────────────────────────────────────

def test_comparison_summary():
    """Test that comparison_summary.py reads and reports compare_runs.py output."""
    tmp = tempfile.mkdtemp(prefix="test_summary_")
    try:
        # Write a summary.tsv and a matched file with a notable delta
        with open(os.path.join(tmp, "summary.tsv"), "w") as fh:
            fh.write("sample\tsnp_old\tsnp_new\tsnp_matched\tsnp_old_only\t"
                     "snp_new_only\tindel_old\tindel_new\tindel_matched\t"
                     "indel_old_only\tindel_new_only\n")
            fh.write("TestSample\t5\t4\t3\t2\t1\t1\t1\t1\t0\t0\n")

        with open(os.path.join(tmp, "TestSample.snp.matched.tsv"), "w") as fh:
            fh.write("Position\tRef\tVarAllele\tOldFreq\tNewFreq\tDeltaFreq\n")
            fh.write("100\tA\tG\t5.0%\t5.1%\t+0.10pp\n")
            fh.write("200\tC\tT\t10.0%\t15.0%\t+5.00pp\n")

        result = subprocess.run(
            [sys.executable, os.path.join(SCRIPT_DIR, "comparison_summary.py"), tmp],
            capture_output=True, text=True)
        assert result.returncode == 0, f"comparison_summary.py failed:\n{result.stderr}"
        assert "TestSample" in result.stdout
        assert "5.00pp" in result.stdout  # notable delta reported
        print("PASS: test_comparison_summary")

    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ── Test: compare_depth_sweep.py ──────────────────────────────────────────────

def test_compare_depth_sweep():
    """Test compare_depth_sweep.py with synthetic depth-sweep consensus files."""
    tmp = tempfile.mkdtemp(prefix="test_depth_sweep_")
    try:
        sample = "TestSample"
        sweep_dir = os.path.join(tmp, "Results", sample, "depth-sweep")
        os.makedirs(sweep_dir)

        # Write consensus files for d=8000 and d=0
        rows_d8000 = [
            {"Position": 10, "Ref": "A", "Var": "G", "VarFreq": "3.5%",
             "Reads1": 7500, "Reads2": 270},
            {"Position": 20, "Ref": "C", "Var": "T", "VarFreq": "2.1%",
             "Reads1": 7800, "Reads2": 167},
            {"Position": 30, "Ref": "G", "Var": "A", "VarFreq": "1.0%",
             "Reads1": 7920, "Reads2": 80},
        ]
        rows_d0 = [
            {"Position": 10, "Ref": "A", "Var": "G", "VarFreq": "3.2%",
             "Reads1": 300000, "Reads2": 9900},
            {"Position": 20, "Ref": "C", "Var": "T", "VarFreq": "2.0%",
             "Reads1": 310000, "Reads2": 6327},
            {"Position": 40, "Ref": "T", "Var": "C", "VarFreq": "1.5%",
             "Reads1": 315000, "Reads2": 4795},
        ]
        _write_varscan_consensus(
            os.path.join(sweep_dir, f"{sample}.d8000.consensus.vcf"), rows_d8000)
        _write_varscan_consensus(
            os.path.join(sweep_dir, f"{sample}.d0.consensus.vcf"), rows_d0)

        out_dir = os.path.join(tmp, "comparison")
        result = subprocess.run(
            [sys.executable, os.path.join(SCRIPT_DIR, "compare_depth_sweep.py"),
             tmp, "--old-depth", "8000", "--new-depth", "0",
             "--output-dir", out_dir],
            capture_output=True, text=True)
        assert result.returncode == 0, f"compare_depth_sweep.py failed:\n{result.stderr}"

        # Verify summary
        summary = _read_tsv(os.path.join(out_dir, "summary.tsv"))
        assert len(summary) == 1
        assert summary[0]["sample"] == sample

        # Verify matched (pos 10, 20 are common)
        matched_path = os.path.join(out_dir, f"{sample}.consensus.matched.tsv")
        assert os.path.isfile(matched_path)
        matched = _read_tsv(matched_path)
        assert len(matched) == 2

        # Verify old_only (pos 30 only in d=8000)
        old_only_path = os.path.join(out_dir, f"{sample}.consensus.old_only.tsv")
        assert os.path.isfile(old_only_path)
        old_only = _read_tsv(old_only_path)
        assert len(old_only) == 1
        assert old_only[0]["Position"] == "30"

        # Verify new_only (pos 40 only in d=0)
        new_only_path = os.path.join(out_dir, f"{sample}.consensus.new_only.tsv")
        assert os.path.isfile(new_only_path)
        new_only = _read_tsv(new_only_path)
        assert len(new_only) == 1
        assert new_only[0]["Position"] == "40"

        print("PASS: test_compare_depth_sweep")

    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ── Test: winpath.py ──────────────────────────────────────────────────────────

def test_winpath():
    """Test Windows-to-WSL path conversion."""
    cases = [
        (r"Z:\2026-05 Cindy\P#5 Amplicon NGS\file.fa",
         "/mnt/z/2026-05 Cindy/P#5 Amplicon NGS/file.fa"),
        (r"C:\Users\test\data.txt", "/mnt/c/Users/test/data.txt"),
    ]
    for winpath, expected in cases:
        result = subprocess.run(
            [sys.executable, os.path.join(SCRIPT_DIR, "winpath.py"), winpath],
            capture_output=True, text=True)
        assert result.returncode == 0, f"winpath.py failed on {winpath}"
        got = result.stdout.strip()
        assert got == expected, f"Expected {expected!r}, got {got!r}"
    print("PASS: test_winpath")


# ── Runner ────────────────────────────────────────────────────────────────────

def main():
    test_compare_runs()
    test_comparison_summary()
    test_compare_depth_sweep()
    test_winpath()
    print("\n=== ALL TESTS PASSED ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
