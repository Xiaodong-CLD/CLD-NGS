"""Unmapped Read Analysis pipeline module for the pipeline Streamlit app.

⚠️  This pipeline integration has not been fully end-to-end tested.

Wraps **unmap-pip** (unmapped read analysis against plasmid + virus databases
using STAR alignment).

Entry point : unmap-pip/scripts/analysis/unmap_analysis.sh
              (preceded by a STAR genome-build step from plasmid .fa files)

Inputs supplied by the user
----------------------------
    BAM directory       Flat directory of input .bam files.
    Plasmid FA dir      Directory of plasmid .fa files. A STAR index is built
                        here under plasmid_star_dir/ before the pipeline runs.
    Virus genome dir    Pre-built STAR virus genome directory.

Output location
---------------
    Results are written to <bam_dir>/../Result/<sample_name>/  by
    unmap_analysis.sh.  The app records <bam_dir>/../Result/ as the output
    directory for job-history purposes.
"""

import pathlib
import shlex
import subprocess

import streamlit as st

from utils import normalize_path, build_output_dir, OUTPUT_ROOT

PIPELINE_ID  = "unmap-pip"
DISPLAY_NAME = "Unmapped Read Analysis (unmap-pip)"

_CONDA_ENV       = "unmap_analysis"
_PIPELINE_DIR    = pathlib.Path(__file__).parents[2] / "unmap-pip"
_PIPELINE_SCRIPT = _PIPELINE_DIR / "scripts" / "analysis" / "unmap_analysis.sh"

# Widget keys for form pre-fill via st.session_state.
_FORM_KEYS = {
    "bam_dir":          "unmap_bam_dir",
    "plasmid_fa_dir":   "unmap_plasmid_fa_dir",
    "virus_genome_dir": "unmap_virus_genome_dir",
    "run_name":         "unmap_run_name",
}


def collect_export_files(output_dir):
    """Return all files under output_dir for export (no intermediate filtering)."""
    out = pathlib.Path(output_dir)
    if not out.exists():
        return []
    return [p for p in out.rglob("*") if p.is_file()]


def _get_conda_base() -> str:
    """Return the conda base directory. Raises ValueError if not found."""
    result = subprocess.run(
        ["conda", "info", "--base"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or not result.stdout.strip():
        raise ValueError(
            "Could not locate conda installation (is conda on PATH?).\n"
            + result.stderr.strip()
        )
    return result.stdout.strip()


def _check_conda_env(conda_env: str) -> None:
    """Raise ValueError if the named conda env is not found."""
    result = subprocess.run(
        ["conda", "env", "list"],
        capture_output=True,
        text=True,
    )
    env_names = [
        line.split()[0]
        for line in result.stdout.splitlines()
        if line and not line.startswith("#")
    ]
    if conda_env not in env_names:
        raise ValueError(
            f"Conda env '{conda_env}' not found.\n"
            f"Create it first: conda env create -f {_PIPELINE_DIR}/envs/unmap_analysis.yml"
        )


def render_form(prefill: dict | None = None) -> dict | None:
    """Render the unmap-pip submission form.

    Returns submitted raw values as a dict, or None before submission.
    """
    if prefill:
        for field, key in _FORM_KEYS.items():
            val = prefill.get(field)
            if val is not None:
                st.session_state[key] = val

    with st.form("unmap_form"):
        st.warning(
            "⚠️ **Not fully tested** — This pipeline integration has not been "
            "end-to-end validated inside the app. Verify outputs manually after the run."
        )

        st.subheader("Input files")

        bam_dir = st.text_input(
            "BAM directory",
            key=_FORM_KEYS["bam_dir"],
            placeholder="/mnt/z/unmap/run1/samples_bam_data",
            help="Directory containing the input .bam files.",
        )
        plasmid_fa_dir = st.text_input(
            "Plasmid FA directory",
            key=_FORM_KEYS["plasmid_fa_dir"],
            placeholder="/mnt/z/unmap/plasmids_data",
            help=(
                "Directory of plasmid .fa files. All .fa files will be concatenated "
                "into plasmid.fa and a STAR genome index will be built under "
                "plasmid_star_dir/ before the analysis runs."
            ),
        )
        virus_genome_dir = st.text_input(
            "Virus genome directory (pre-built STAR index)",
            key=_FORM_KEYS["virus_genome_dir"],
            placeholder="/mnt/z/unmap/viral_db/unmasked_virus_dir_v3",
            help="Pre-built STAR virus genome directory.",
        )

        st.subheader("Run name")
        run_name = st.text_input(
            "Run name",
            key=_FORM_KEYS["run_name"],
            placeholder="e.g. unmap_run1",
            help=(
                "Short label for tracking this run in the app (no slashes). "
                "Results are written directly to <bam_dir>/../Result/ by the pipeline."
            ),
        )

        submitted = st.form_submit_button("Submit")

    if not submitted:
        return None

    return {
        "bam_dir":          bam_dir,
        "plasmid_fa_dir":   plasmid_fa_dir,
        "virus_genome_dir": virus_genome_dir,
        "run_name":         run_name,
    }


def parse_command(cmd_str: str) -> dict:
    """Parse a previous command string into a prefill dict.

    unmap commands are multi-step bash -c '...' strings; structured parsing
    is not feasible.  Returns empty strings so the user fills in the form.
    """
    return {
        "bam_dir":          "",
        "plasmid_fa_dir":   "",
        "virus_genome_dir": "",
        "run_name":         "",
    }


def validate_and_build_command(
    raw_values: dict,
) -> tuple[dict, list[str], str]:
    """Validate inputs and build the pipeline command.

    Dependency checks (conda env existence) are run synchronously here so
    the user sees a clear error before a job is submitted.

    Returns:
        normalized_config : dict   Values saved with the job history entry.
        command           : list   Pipeline command as an argument list.
        working_directory : str    Directory from which the command runs.

    Raises ValueError with a user-facing message on invalid input.
    """
    for field, label in [
        ("bam_dir",          "BAM directory"),
        ("plasmid_fa_dir",   "Plasmid FA directory"),
        ("virus_genome_dir", "Virus genome directory"),
        ("run_name",         "Run name"),
    ]:
        if not raw_values.get(field, "").strip():
            raise ValueError(f"{label} is required.")

    run_name = raw_values["run_name"].strip()
    if any(c in run_name for c in ("/", "\\")):
        raise ValueError("Run name must not contain path separators (/ or \\).")

    bam_dir          = normalize_path(raw_values["bam_dir"])
    plasmid_fa_dir   = normalize_path(raw_values["plasmid_fa_dir"])
    virus_genome_dir = normalize_path(raw_values["virus_genome_dir"])

    if not bam_dir.is_dir():
        raise ValueError(f"BAM directory not found: {bam_dir}")
    if not plasmid_fa_dir.is_dir():
        raise ValueError(f"Plasmid FA directory not found: {plasmid_fa_dir}")
    if not virus_genome_dir.is_dir():
        raise ValueError(f"Virus genome directory not found: {virus_genome_dir}")
    if not _PIPELINE_SCRIPT.is_file():
        raise ValueError(
            f"unmap_analysis.sh not found: {_PIPELINE_SCRIPT}\n"
            "Is the unmap-pip directory present?"
        )

    # Dependency checks — fail fast before submitting the job.
    conda_base = _get_conda_base()
    conda_sh   = pathlib.Path(conda_base) / "etc" / "profile.d" / "conda.sh"
    if not conda_sh.is_file():
        raise ValueError(
            f"conda.sh not found at expected location: {conda_sh}\n"
            "Check your conda installation."
        )
    _check_conda_env(_CONDA_ENV)

    plasmid_genome_dir = plasmid_fa_dir / "plasmid_star_dir"
    # unmap_analysis.sh writes to <bam_dir>/../Result/<sample>/
    output_dir = bam_dir.parent / "Result"

    # Build multi-step command that mirrors job.script:
    #   1. source conda & activate env
    #   2. build STAR plasmid index
    #   3. run unmap_analysis.sh
    inner_cmd = " && ".join([
        f"source {shlex.quote(str(conda_sh))}",
        f"conda activate {shlex.quote(_CONDA_ENV)}",
        f"mkdir -p {shlex.quote(str(plasmid_genome_dir))}",
        (
            f"cat {shlex.quote(str(plasmid_fa_dir))}/*.fa "
            f"> {shlex.quote(str(plasmid_fa_dir / 'plasmid.fa'))}"
        ),
        (
            f"STAR --runMode genomeGenerate "
            f"--genomeDir {shlex.quote(str(plasmid_genome_dir))} "
            f"--genomeFastaFiles {shlex.quote(str(plasmid_fa_dir / 'plasmid.fa'))} "
            f"--runThreadN 24"
        ),
        (
            f"bash {shlex.quote(str(_PIPELINE_SCRIPT))} "
            f"{shlex.quote(str(bam_dir))} "
            f"{shlex.quote(str(plasmid_genome_dir))} "
            f"{shlex.quote(str(virus_genome_dir))}"
        ),
    ])
    command = ["bash", "-c", inner_cmd]

    normalized_config = {
        "bam_dir":          str(bam_dir),
        "plasmid_fa_dir":   str(plasmid_fa_dir),
        "virus_genome_dir": str(virus_genome_dir),
        "run_name":         run_name,
        "output_dir":       str(output_dir),
    }

    return normalized_config, command, str(_PIPELINE_DIR)
