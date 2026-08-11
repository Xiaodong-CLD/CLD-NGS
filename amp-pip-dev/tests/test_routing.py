"""Smoke test: verify that pipeline_combined.py routes intermediates to
tmp/<sample>/ and copies final outputs to Results/<sample>/.

Run from amp-pip-dev/:

    cd amp-pip-dev
    conda run -n variant-calling pytest tests/test_routing.py -v

What is tested:
  - FASTQ symlinks are placed in tmp/<sample>/, not Results/<sample>/.
  - tmp/<sample>/ and Results/<sample>/ are both created.
  - fastp QC reports, coverage plots, varscan SNP/indel VCFs, and the
    pindel Summary VCF are all copied to Results/<sample>/.
  - Run-level Results/snp.txt and Results/indel.txt are written.
  - No BAM / mpileup files appear under Results/.

The test mocks subprocess.call so bwa / samtools / pindel / fastp / varscan
never run. The files they would have written are pre-seeded in tmp/<sample>/
by the mock's side-effect, mimicking a real pipeline run.
"""

import pathlib
import shlex
import sys
import unittest.mock as mock

import openpyxl
import pytest

_AMP_PIP_DEV = pathlib.Path(__file__).parents[1].resolve()
_SCRIPT = _AMP_PIP_DEV / "pipeline_combined.py"

# The sample sheet uses naming convention "NAME_Lxxx_…" → naming() returns "NAME".
SAMPLE_ID   = "S1"          # xlsx row index (used as SampleID in report)
FQ_FWD      = "S1_L001_R1_001.fastq.gz"
FQ_REV      = "S1_L001_R2_001.fastq.gz"
REFERENCE   = "goi.fa"
SAMPLE_NAME = "S1"          # naming("S1_L001_R1_001.fastq.gz") → "S1"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def run_dir(tmp_path):
    """Create a minimal amp-pip-dev input directory layout."""
    d = tmp_path / "test_run"
    (d / "Fastq").mkdir(parents=True)
    (d / "reference").mkdir()
    (d / "sample_xlsx").mkdir()

    # Minimal FASTQ pair (content irrelevant — pipeline is mocked).
    (d / "Fastq" / FQ_FWD).write_bytes(b"@r\nACGT\n+\nIIII\n")
    (d / "Fastq" / FQ_REV).write_bytes(b"@r\nACGT\n+\nIIII\n")

    # Reference FASTA.
    (d / "reference" / REFERENCE).write_text(">seq\nACGT\n")

    # Sample spreadsheet: one sample, one clone.
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["", "fastq (forward)", "fastq (reverse)",
               "GOI reference file", "Clone ID (optional)"])
    ws.append([SAMPLE_ID, FQ_FWD, FQ_REV, REFERENCE, "clone1"])
    wb.save(d / "sample_xlsx" / "samples.xlsx")
    return d


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _seed_tmp_outputs(tmp_sample: pathlib.Path) -> None:
    """Create the files that a real pipeline.py + amino-anno run would produce
    under tmp/<sample>/.  Called by the subprocess.call mock when pipeline.py
    is invoked.
    """
    tmp_sample.mkdir(parents=True, exist_ok=True)

    # fastp QC reports (written alongside the input FASTQ symlinks).
    (tmp_sample / f"{SAMPLE_NAME}.fastp.html").write_text("<html/>")
    (tmp_sample / f"{SAMPLE_NAME}.fastp.json").write_text("{}")

    # Coverage-plot output directory.
    cov = tmp_sample / "coverage-plot-output"
    cov.mkdir(exist_ok=True)
    (cov / f"{SAMPLE_NAME}_cov.html").write_text("<html/>")

    # VarScan outputs.
    vs = tmp_sample / "varscan-output"
    vs.mkdir(exist_ok=True)
    (vs / f"{SAMPLE_NAME}.snp.vcf").write_text("# snp vcf\n")
    (vs / f"{SAMPLE_NAME}.indel.vcf").write_text("# indel vcf\n")
    # Amino-annotated SNP file: header-only TSV so pd.read_csv succeeds
    # (all processing on an empty DataFrame is safe).
    aa_cols = [
        "Chrom", "Position", "Cons", "Reads1", "Reads2", "VarFreq",
        "Strands1", "Strands2", "Qual1", "Qual2",
        "MapQual1", "MapQual2", "Pvalue", "MapPvalue",
        "Reads1Plus", "Reads1Minus", "Reads2Plus", "Reads2Minus", "VarAllele",
    ]
    (vs / f"{SAMPLE_NAME}.snp.amino_annotated.vcf").write_text(
        "\t".join(aa_cols) + "\n"
    )

    # Pindel Summary VCF: VCF-format header only, no data rows.
    p = tmp_sample / "pindel-output"
    p.mkdir(exist_ok=True)
    (p / "pindel-output.Summary.vcf").write_text("##fileformat=VCFv4.1\n")

    # LoFreq indel VCF: VCF-format header only, no data rows.
    lf = tmp_sample / "lofreq-output"
    lf.mkdir(exist_ok=True)
    (lf / f"{SAMPLE_NAME}.lofreq.indel.vcf").write_text("##fileformat=VCFv4.1\n")


def _make_subprocess_side_effect():
    """Return a side-effect for subprocess.call that:
    - seeds tmp/<sample>/ when pipeline.py is invoked,
    - is a no-op for all other subprocess.call invocations.
    """
    def _side_effect(cmd, *args, **kwargs):
        parts = cmd if isinstance(cmd, list) else shlex.split(str(cmd))
        if any("pipeline.py" in str(p) for p in parts):
            # Locate the -f1 argument to find the tmp/<sample>/ directory.
            for i, part in enumerate(parts):
                if str(part) in ("-f1", "--input-file-1") and i + 1 < len(parts):
                    tmp_sample = pathlib.Path(str(parts[i + 1])).parent
                    _seed_tmp_outputs(tmp_sample)
                    break
        return 0
    return _side_effect


# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------

def test_tmp_routing_and_results_copy(run_dir):
    """Intermediates land in tmp/<sample>/ and final outputs are copied to
    Results/<sample>/ by pipeline_combined.py."""
    import runpy

    with mock.patch("subprocess.call",
                    side_effect=_make_subprocess_side_effect()), \
         mock.patch("sys.argv",
                    ["pipeline_combined.py", str(run_dir)]):
        runpy.run_path(str(_SCRIPT), run_name="__main__")

    tmp_s = run_dir / "tmp" / SAMPLE_NAME
    res_s = run_dir / "Results" / SAMPLE_NAME

    # ── Intermediates land in tmp/<sample>/ ──────────────────────────────────
    assert tmp_s.is_dir(), "tmp/<sample>/ was not created"
    assert (tmp_s / f"{SAMPLE_NAME}.R1.fastq.gz").is_symlink(), \
        "FASTQ symlink was not placed in tmp/<sample>/"

    # ── Final outputs copied to Results/<sample>/ ────────────────────────────
    assert res_s.is_dir(), "Results/<sample>/ was not created"
    assert (res_s / f"{SAMPLE_NAME}.fastp.html").is_file(), \
        "fastp HTML missing from Results/<sample>/"
    assert (res_s / "coverage-plot-output" / f"{SAMPLE_NAME}_cov.html").is_file(), \
        "coverage plot missing from Results/<sample>/coverage-plot-output/"
    assert (res_s / "varscan-output" / f"{SAMPLE_NAME}.snp.amino_annotated.vcf").is_file(), \
        "annotated SNP VCF missing from Results/<sample>/varscan-output/"
    assert (res_s / "varscan-output" / f"{SAMPLE_NAME}.indel.vcf").is_file(), \
        "indel VCF missing from Results/<sample>/varscan-output/"
    assert (res_s / "pindel-output" / "pindel-output.Summary.vcf").is_file(), \
        "pindel Summary VCF missing from Results/<sample>/pindel-output/"
    assert (res_s / "lofreq-output" / f"{SAMPLE_NAME}.lofreq.indel.vcf").is_file(), \
        "lofreq indel VCF missing from Results/<sample>/lofreq-output/"

    # ── Run-level combined reports exist in Results/ ─────────────────────────
    assert (run_dir / "Results" / "snp.txt").is_file(), "snp.txt missing"
    assert (run_dir / "Results" / "indel.txt").is_file(), "indel.txt missing"

    # ── Intermediates did NOT leak into Results/ ─────────────────────────────
    assert not any(res_s.rglob("*.bam")), "BAM file leaked into Results/"
    assert not any(res_s.rglob("*.mpileup")), "mpileup leaked into Results/"
