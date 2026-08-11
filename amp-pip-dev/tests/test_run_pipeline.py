"""End-to-end smoke test for run_pipeline.py.

Verifies that the compatibility layer (run_pipeline.py →
prepare_pipeline_input.sh → pipeline_combined.py) produces the same
directory structure and output files as calling pipeline_combined.py
directly with the pre-built input layout.

Run from amp-pip-dev/:

    cd amp-pip-dev
    conda run -n variant-calling pytest tests/test_run_pipeline.py -v

What is tested:
  - run_pipeline.py calls prepare_pipeline_input.sh with the correct four
    arguments (fastq_dir, sample_sheet, reference, output_dir).
  - The expected symlink structure (Fastq/, sample_xlsx/, reference/) is
    created inside the output directory.
  - pipeline_combined.py is called with the output directory as input_dir.
  - The final Results/ layout matches what test_routing.py expects from
    pipeline_combined.py directly (FASTQ symlinks in tmp/, VCFs and reports
    copied to Results/, snp.txt and indel.txt written at run level).

subprocess.run is mocked for both steps so no real tools run:
  - The bash-script call is replaced by a Python helper that creates the
    same symlinks the script would create.
  - The pipeline_combined.py call is replaced by running the script in-
    process via runpy with subprocess.call mocked to seed tool outputs.
"""

import os
import pathlib
import runpy
import shlex
import sys
import unittest.mock as mock

import openpyxl
import pytest

_AMP_PIP_DEV = pathlib.Path(__file__).parents[1].resolve()
_RUN_PIPELINE = _AMP_PIP_DEV / "run_pipeline.py"
_PIPELINE_COMBINED = _AMP_PIP_DEV / "pipeline_combined.py"
_PREPARE_SH = _AMP_PIP_DEV / "tools" / "prepare_pipeline_input.sh"

SAMPLE_ID   = "S1"
FQ_FWD      = "S1_L001_R1_001.fastq.gz"
FQ_REV      = "S1_L001_R2_001.fastq.gz"
REFERENCE   = "goi.fa"
SAMPLE_NAME = "S1"   # naming("S1_L001_R1_001.fastq.gz") → "S1"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def user_inputs(tmp_path):
    """Create separate source directories matching user-supplied paths."""
    fastq_dir = tmp_path / "fastq_source"
    fastq_dir.mkdir()
    (fastq_dir / FQ_FWD).write_bytes(b"@r\nACGT\n+\nIIII\n")
    (fastq_dir / FQ_REV).write_bytes(b"@r\nACGT\n+\nIIII\n")

    ref_dir = tmp_path / "ref_source"
    ref_dir.mkdir()
    (ref_dir / REFERENCE).write_text(">seq\nACGT\n")

    # Sample sheet in its own directory (user provides the file path).
    sheet_dir = tmp_path / "sheets"
    sheet_dir.mkdir()
    sheet_path = sheet_dir / "samples.xlsx"
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["", "fastq (forward)", "fastq (reverse)",
               "GOI reference file", "Clone ID (optional)"])
    ws.append([SAMPLE_ID, FQ_FWD, FQ_REV, REFERENCE, "clone1"])
    wb.save(sheet_path)

    output_dir = tmp_path / "output"

    return {
        "fastq_dir": fastq_dir,
        "sample_sheet": sheet_path,
        "reference": ref_dir / REFERENCE,
        "output_dir": output_dir,
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _simulate_prepare_sh(fastq_dir, sample_sheet, reference, target_dir):
    """Pure-Python equivalent of prepare_pipeline_input.sh.

    Creates the same Fastq/, sample_xlsx/, reference/ symlink layout that
    the bash script would create, without requiring bash or realpath.
    """
    fastq_dir   = pathlib.Path(fastq_dir).resolve()
    sample_sheet = pathlib.Path(sample_sheet).resolve()
    reference   = pathlib.Path(reference).resolve()
    target_dir  = pathlib.Path(target_dir)

    (target_dir / "Fastq").mkdir(parents=True, exist_ok=True)
    (target_dir / "sample_xlsx").mkdir(exist_ok=True)
    (target_dir / "reference").mkdir(exist_ok=True)

    for src in sorted(fastq_dir.iterdir()):
        if src.suffix in {".fastq", ".gz", ".fq"} or src.name.endswith(
            (".fastq.gz", ".fq.gz")
        ):
            dst = target_dir / "Fastq" / src.name
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            os.symlink(src, dst)

    sheet_dst = target_dir / "sample_xlsx" / sample_sheet.name
    if sheet_dst.is_symlink() or sheet_dst.exists():
        sheet_dst.unlink()
    os.symlink(sample_sheet, sheet_dst)

    ref_dst = target_dir / "reference" / reference.name
    if ref_dst.is_symlink() or ref_dst.exists():
        ref_dst.unlink()
    os.symlink(reference, ref_dst)


def _seed_tmp_outputs(tmp_sample: pathlib.Path) -> None:
    """Seed files that pipeline.py + amino-anno would write under tmp/<sample>/."""
    tmp_sample.mkdir(parents=True, exist_ok=True)

    (tmp_sample / f"{SAMPLE_NAME}.fastp.html").write_text("<html/>")
    (tmp_sample / f"{SAMPLE_NAME}.fastp.json").write_text("{}")

    cov = tmp_sample / "coverage-plot-output"
    cov.mkdir(exist_ok=True)
    (cov / f"{SAMPLE_NAME}_cov.html").write_text("<html/>")

    vs = tmp_sample / "varscan-output"
    vs.mkdir(exist_ok=True)
    (vs / f"{SAMPLE_NAME}.snp.vcf").write_text("# snp vcf\n")
    (vs / f"{SAMPLE_NAME}.indel.vcf").write_text("# indel vcf\n")

    aa_cols = [
        "Chrom", "Position", "Cons", "Reads1", "Reads2", "VarFreq",
        "Strands1", "Strands2", "Qual1", "Qual2",
        "MapQual1", "MapQual2", "Pvalue", "MapPvalue",
        "Reads1Plus", "Reads1Minus", "Reads2Plus", "Reads2Minus", "VarAllele",
    ]
    (vs / f"{SAMPLE_NAME}.snp.amino_annotated.vcf").write_text(
        "\t".join(aa_cols) + "\n"
    )

    p = tmp_sample / "pindel-output"
    p.mkdir(exist_ok=True)
    (p / "pindel-output.Summary.vcf").write_text("##fileformat=VCFv4.1\n")

    lf = tmp_sample / "lofreq-output"
    lf.mkdir(exist_ok=True)
    (lf / f"{SAMPLE_NAME}.lofreq.indel.vcf").write_text("##fileformat=VCFv4.1\n")


def _make_pipeline_call_side_effect():
    """subprocess.call mock for pipeline_combined.py's tool invocations."""
    def _side_effect(cmd, *args, **kwargs):
        parts = cmd if isinstance(cmd, list) else shlex.split(str(cmd))
        if any("pipeline.py" in str(p) for p in parts):
            for i, part in enumerate(parts):
                if str(part) in ("-f1", "--input-file-1") and i + 1 < len(parts):
                    tmp_sample = pathlib.Path(str(parts[i + 1])).parent
                    _seed_tmp_outputs(tmp_sample)
                    break
        return 0
    return _side_effect


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_run_pipeline_calls_prepare_sh_with_correct_args(user_inputs):
    """run_pipeline.py must invoke prepare_pipeline_input.sh with the four
    user-supplied paths in the correct order."""
    captured_calls = []

    def _mock_run(cmd, check=False, **kw):
        captured_calls.append(list(cmd))
        # Short-circuit both calls; we only check arguments here.
        return mock.MagicMock(returncode=0)

    with mock.patch("subprocess.run", side_effect=_mock_run), \
         mock.patch("sys.argv", [
             "run_pipeline.py",
             "--sample-sheet", str(user_inputs["sample_sheet"]),
             "--fastq-dir",    str(user_inputs["fastq_dir"]),
             "--reference",    str(user_inputs["reference"]),
             "--output-dir",   str(user_inputs["output_dir"]),
         ]):
        runpy.run_path(str(_RUN_PIPELINE), run_name="__main__")

    assert len(captured_calls) == 2, "Expected exactly two subprocess.run calls"

    # First call: bash prepare_pipeline_input.sh <fastq_dir> <sheet> <ref> <output>
    prepare_call = captured_calls[0]
    assert prepare_call[0] == "bash"
    assert str(_PREPARE_SH) in prepare_call[1]
    assert str(user_inputs["fastq_dir"])   in prepare_call[2]
    assert str(user_inputs["sample_sheet"]) in prepare_call[3]
    assert str(user_inputs["reference"])   in prepare_call[4]
    assert str(user_inputs["output_dir"])  in prepare_call[5]

    # Second call: python pipeline_combined.py <output_dir>
    combined_call = captured_calls[1]
    assert str(_PIPELINE_COMBINED) in combined_call[1]
    assert str(user_inputs["output_dir"]) in combined_call[2]


def test_run_pipeline_produces_same_results_as_direct_call(user_inputs):
    """The compatibility layer (run_pipeline.py) produces Results/snp.txt,
    Results/indel.txt, and per-sample outputs identical in structure to
    calling pipeline_combined.py directly with the pre-built layout."""
    output_dir = user_inputs["output_dir"]

    def _mock_run(cmd, check=False, **kw):
        cmd_list = list(cmd)
        str_cmd = " ".join(str(c) for c in cmd_list)
        if "prepare_pipeline_input" in str_cmd:
            # Simulate the bash script creating symlinks.
            _simulate_prepare_sh(
                user_inputs["fastq_dir"],
                user_inputs["sample_sheet"],
                user_inputs["reference"],
                output_dir,
            )
        elif "pipeline_combined" in str_cmd:
            # Run pipeline_combined.py in-process with tool calls mocked.
            with mock.patch("subprocess.call",
                            side_effect=_make_pipeline_call_side_effect()), \
                 mock.patch("sys.argv",
                            ["pipeline_combined.py", str(output_dir)]):
                runpy.run_path(str(_PIPELINE_COMBINED), run_name="__main__")
        return mock.MagicMock(returncode=0)

    with mock.patch("subprocess.run", side_effect=_mock_run), \
         mock.patch("sys.argv", [
             "run_pipeline.py",
             "--sample-sheet", str(user_inputs["sample_sheet"]),
             "--fastq-dir",    str(user_inputs["fastq_dir"]),
             "--reference",    str(user_inputs["reference"]),
             "--output-dir",   str(output_dir),
         ]):
        runpy.run_path(str(_RUN_PIPELINE), run_name="__main__")

    # ── Symlink structure created inside output_dir ──────────────────────────
    assert (output_dir / "Fastq" / FQ_FWD).is_symlink(), \
        "FASTQ symlink missing from output_dir/Fastq/"
    assert (output_dir / "sample_xlsx" / user_inputs["sample_sheet"].name).is_symlink(), \
        "sample sheet symlink missing from output_dir/sample_xlsx/"
    assert (output_dir / "reference" / REFERENCE).is_symlink(), \
        "reference symlink missing from output_dir/reference/"

    # ── FASTQ renames land in tmp/, not Results/ ─────────────────────────────
    tmp_s = output_dir / "tmp" / SAMPLE_NAME
    assert tmp_s.is_dir(), "tmp/<sample>/ was not created"
    assert (tmp_s / f"{SAMPLE_NAME}.R1.fastq.gz").is_symlink(), \
        "renamed FASTQ symlink missing from tmp/<sample>/"

    # ── Final outputs in Results/<sample>/ ───────────────────────────────────
    res_s = output_dir / "Results" / SAMPLE_NAME
    assert res_s.is_dir(), "Results/<sample>/ was not created"
    assert (res_s / f"{SAMPLE_NAME}.fastp.html").is_file()
    assert (res_s / "coverage-plot-output" / f"{SAMPLE_NAME}_cov.html").is_file()
    assert (res_s / "varscan-output" / f"{SAMPLE_NAME}.snp.amino_annotated.vcf").is_file()
    assert (res_s / "varscan-output" / f"{SAMPLE_NAME}.indel.vcf").is_file()
    assert (res_s / "pindel-output" / "pindel-output.Summary.vcf").is_file()
    assert (res_s / "lofreq-output" / f"{SAMPLE_NAME}.lofreq.indel.vcf").is_file()

    # ── Run-level reports ────────────────────────────────────────────────────
    assert (output_dir / "Results" / "snp.txt").is_file(), "snp.txt missing"
    assert (output_dir / "Results" / "indel.txt").is_file(), "indel.txt missing"

    # ── Intermediates did not leak into Results/ ─────────────────────────────
    assert not any(res_s.rglob("*.bam")),     "BAM leaked into Results/"
    assert not any(res_s.rglob("*.mpileup")), "mpileup leaked into Results/"
