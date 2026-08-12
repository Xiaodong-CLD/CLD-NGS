"""Tests for the desktop-export helpers in utils.py.

Run from the pipeline_app directory:

    cd pipeline_app
    python -m pytest tests/test_export.py -v

All file operations use pytest's tmp_path — nothing is written to the real
Windows desktop.
"""

import os
import zipfile

import pytest

import utils


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_output_dir(root):
    """Create a synthetic pipeline output dir mixing final results and junk."""
    out = root / "260803_run"
    # Final combined reports.
    fr = out / "final_report"
    (fr / "coverage_plots").mkdir(parents=True)
    (fr / "all_snv.tsv").write_text("snv\n")
    (fr / "all_indel.tsv").write_text("indel\n")
    (fr / "all_consensus.fa").write_text(">c\nACGT\n")
    (fr / "coverage_plots" / "sampleA_coverage_plot.jpeg").write_bytes(b"\xff\xd8jpeg")

    # Per-sample final results (Results/), including a large IGV BAM to exclude.
    res = out / "Results" / "sampleA"
    (res / "snv").mkdir(parents=True)
    (res / "consensus").mkdir(parents=True)
    (res / "uninformative_analysis" / "igv").mkdir(parents=True)
    (res / "snv" / "varscan_snp_report.tsv").write_text("report\n")
    (res / "consensus" / "sampleA.consensus.fa").write_text(">s\nACGT\n")
    (res / "uninformative_analysis" / "igv" / "sampleA_bamout.bam").write_bytes(b"BAMDATA")
    (res / "uninformative_analysis" / "igv" / "sampleA_bamout.bai").write_bytes(b"BAIDATA")

    # Logs (kept).
    (out / "logs").mkdir()
    (out / "logs" / "pipeline.log").write_text("done\n")

    # Top-level docs/config (kept).
    (out / "README.md").write_text("# results\n")
    (out / "run_params.json").write_text("{}\n")

    # tmp/ — all intermediate, must be excluded entirely.
    tmp = out / "tmp" / "sampleA"
    (tmp).mkdir(parents=True)
    (tmp / "sampleA.sam").write_text("sam\n")
    (tmp / "sampleA.sorted.bam").write_bytes(b"BAM")
    (tmp / "sampleA.R1.trimmed.fastq.gz").write_bytes(b"FASTQ")
    (tmp / "sampleA.mpileup").write_text("pileup\n")
    ref = out / "tmp" / "reference"
    ref.mkdir(parents=True)
    (ref / "amplicon_reference.fa").write_text(">r\nACGT\n")
    (ref / "amplicon_reference.fa.bwt").write_bytes(b"bwt")

    # A stray large intermediate at top level (excluded by suffix).
    (out / "big.bam").write_bytes(b"BAM")
    return out


@pytest.fixture
def output_dir(tmp_path):
    return _make_output_dir(tmp_path)


@pytest.fixture
def desktop(tmp_path):
    d = tmp_path / "Desktop"
    d.mkdir()
    return d


# ---------------------------------------------------------------------------
# collect_final_files — inclusion / exclusion
# ---------------------------------------------------------------------------

def _rel_names(output_dir):
    return {str(rel) for _abs, rel in utils.collect_final_files(output_dir)}


def test_includes_required_final_files(output_dir):
    names = _rel_names(output_dir)
    assert os.path.join("final_report", "all_snv.tsv") in names
    assert os.path.join("final_report", "all_indel.tsv") in names
    assert os.path.join("final_report", "all_consensus.fa") in names
    assert os.path.join("final_report", "coverage_plots", "sampleA_coverage_plot.jpeg") in names
    assert os.path.join("Results", "sampleA", "snv", "varscan_snp_report.tsv") in names
    assert os.path.join("Results", "sampleA", "consensus", "sampleA.consensus.fa") in names
    assert os.path.join("logs", "pipeline.log") in names
    assert "README.md" in names
    assert "run_params.json" in names


def test_excludes_large_and_intermediate_files(output_dir):
    names = _rel_names(output_dir)
    # tmp/ tree fully excluded.
    assert not any(n.startswith("tmp" + os.sep) for n in names)
    # BAM/SAM/FASTQ/pileup and reference indexes excluded everywhere.
    assert not any(n.endswith((".bam", ".bai", ".sam", ".mpileup", ".bwt")) for n in names)
    assert not any(n.endswith(".fastq.gz") for n in names)
    assert "big.bam" not in names
    # IGV bam inside a kept Results dir is still excluded by suffix.
    assert os.path.join(
        "Results", "sampleA", "uninformative_analysis", "igv", "sampleA_bamout.bam"
    ) not in names


def test_empty_output_dir_yields_no_files(tmp_path):
    empty = tmp_path / "empty_run"
    empty.mkdir()
    assert utils.collect_final_files(empty) == []


# ---------------------------------------------------------------------------
# sanitize_output_name
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("my_run", "my_run"),
    ("F8KO pilot", "F8KO pilot"),
    ('bad:name/with\\chars?*', "bad_name_with_chars__"),
    ("trailing.dot.", "trailing.dot"),
    ("  spaced  ", "spaced"),
    ("results.zip", "results"),
    ("", "results"),
    ("CON", "_CON"),
])
def test_sanitize_output_name(raw, expected):
    assert utils.sanitize_output_name(raw) == expected


# ---------------------------------------------------------------------------
# export_results_zip
# ---------------------------------------------------------------------------

def test_creates_pipeline_subfolder_and_zip(output_dir, desktop):
    info = utils.export_results_zip(
        output_dir=str(output_dir),
        pipeline_id="ko-pip",
        output_name="my_run",
        desktop_wsl_path=str(desktop),
    )
    zip_dir = desktop / utils.EXPORT_FOLDER_NAME / "ko-pip"
    assert zip_dir.is_dir()
    zip_path = zip_dir / "my_run.zip"
    assert zip_path.is_file()
    assert info["zip_path"] == str(zip_path)
    assert info["size_bytes"] > 0
    assert info["file_count"] > 0


def test_zip_contains_only_final_files(output_dir, desktop):
    utils.export_results_zip(
        output_dir=str(output_dir),
        pipeline_id="ko-pip",
        output_name="my_run",
        desktop_wsl_path=str(desktop),
    )
    zip_path = desktop / utils.EXPORT_FOLDER_NAME / "ko-pip" / "my_run.zip"
    with zipfile.ZipFile(zip_path) as zf:
        arcnames = zf.namelist()
    # Everything is nested under the <output_name>/ folder.
    assert all(a.startswith("my_run/") for a in arcnames)
    joined = "\n".join(arcnames)
    assert "final_report/all_snv.tsv" in joined
    assert "logs/pipeline.log" in joined
    # No excluded artifacts.
    assert not any(a.endswith((".bam", ".sam", ".fastq.gz", ".mpileup", ".bwt")) for a in arcnames)
    assert "big.bam" not in joined


def test_no_silent_overwrite(output_dir, desktop):
    kwargs = dict(
        output_dir=str(output_dir),
        pipeline_id="ko-pip",
        output_name="my_run",
        desktop_wsl_path=str(desktop),
    )
    utils.export_results_zip(**kwargs)
    # Second export without overwrite must refuse.
    with pytest.raises(utils.ExportZipExists):
        utils.export_results_zip(**kwargs)
    # With overwrite=True it succeeds.
    info = utils.export_results_zip(overwrite=True, **kwargs)
    assert info["file_count"] > 0


def test_output_name_is_sanitized_in_zip_filename(output_dir, desktop):
    info = utils.export_results_zip(
        output_dir=str(output_dir),
        pipeline_id="ko-pip",
        output_name="bad:name?",
        desktop_wsl_path=str(desktop),
    )
    assert info["zip_path"].endswith("bad_name_.zip")


def test_missing_final_files_raises(tmp_path, desktop):
    empty = tmp_path / "empty_run"
    empty.mkdir()
    with pytest.raises(utils.ExportError):
        utils.export_results_zip(
            output_dir=str(empty),
            pipeline_id="ko-pip",
            output_name="run",
            desktop_wsl_path=str(desktop),
        )


def test_missing_output_dir_raises(tmp_path, desktop):
    with pytest.raises(utils.ExportError):
        utils.export_results_zip(
            output_dir=str(tmp_path / "does_not_exist"),
            pipeline_id="ko-pip",
            output_name="run",
            desktop_wsl_path=str(desktop),
        )


def test_unwritable_desktop_raises(output_dir, tmp_path):
    ro_desktop = tmp_path / "ReadOnlyDesktop"
    ro_desktop.mkdir()
    os.chmod(ro_desktop, 0o500)  # read + execute, no write
    try:
        with pytest.raises(utils.ExportError):
            utils.export_results_zip(
                output_dir=str(output_dir),
                pipeline_id="ko-pip",
                output_name="run",
                desktop_wsl_path=str(ro_desktop),
            )
    finally:
        os.chmod(ro_desktop, 0o700)  # restore so tmp cleanup can remove it
