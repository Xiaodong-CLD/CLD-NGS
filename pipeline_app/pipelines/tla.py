"""TLA pipeline module for the pipeline Streamlit app.

Entry point : tla-pip/pipeline.py
Working dir : tla-pip/

Required CLI args:
    -s / --sample-sheet    Sample sheet (.tsv or .xlsx)
    -f / --fastq-dir       Directory containing paired FASTQ files
    -v / --vector-fasta    Vector/GOI FASTA file
    -o / --output-dir      Output directory

Optional CLI args:
    -r / --references-dir  Local host genome + mixed-ref cache (default: auto)
    -a / --accession       NCBI host genome accession (default: GCF_003668045.3)
    --force-download       Re-download host genome even if cached
"""

import pathlib
import subprocess

import streamlit as st

from utils import normalize_path, build_output_dir, OUTPUT_ROOT

PIPELINE_ID = "tla-pip"
DISPLAY_NAME = "TLA Pipeline"

_DEFAULT_ACCESSION = "GCF_003668045.3"
_PIPELINE_DIR = pathlib.Path(__file__).parents[2] / "tla-pip"
_PIPELINE_SCRIPT = _PIPELINE_DIR / "pipeline.py"
_CONDA_ENV = "tla_pip"

# Widget keys for form pre-fill via st.session_state.
_FORM_KEYS = {
    "sample_sheet":   "tla_sample_sheet",
    "fastq_dir":      "tla_fastq_dir",
    "vector_fasta":   "tla_vector_fasta",
    "run_name":       "tla_run_name",
    "references_dir": "tla_references_dir",
    "accession":      "tla_accession",
}


def _tla_python() -> str:
    """Return the Python executable from the tla_pip conda environment.

    Raises ValueError with a clear message if the environment is not found.
    """
    result = subprocess.run(
        ["conda", "run", "-n", _CONDA_ENV, "python", "-c",
         "import sys; print(sys.executable)"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or not result.stdout.strip():
        raise ValueError(
            f"Could not locate Python in conda env '{_CONDA_ENV}'.\n"
            f"Create it first: conda env create -f tla-pip/envs/tla_pip.yml\n"
            f"{result.stderr.strip()}"
        )
    return result.stdout.strip()


def render_form(prefill: dict | None = None) -> dict | None:
    """Render the TLA submission form.

    If *prefill* is provided, widgets are pre-populated so the user only needs
    to update the output directory before submitting.

    Returns submitted raw values as a dict, or None before submission.
    """
    if prefill:
        for field, key in _FORM_KEYS.items():
            val = prefill.get(field)
            if val is not None:
                st.session_state[key] = val

    with st.form("tla_form"):
        st.subheader("Required inputs")

        sample_sheet = st.text_input(
            "Sample sheet (.xlsx)",
            key=_FORM_KEYS["sample_sheet"],
            placeholder=r"D:\TLA\run1\samples.xlsx  or  /mnt/d/TLA/run1/samples.xlsx",
            help="Path to the sample sheet. Windows and WSL paths are both accepted.",
        )
        fastq_dir = st.text_input(
            "FASTQ directory",
            key=_FORM_KEYS["fastq_dir"],
            placeholder=r"D:\TLA\run1\Fastq  or  /mnt/d/TLA/run1/Fastq",
            help="Directory containing the paired FASTQ files listed in the sample sheet.",
        )
        vector_fasta = st.text_input(
            "Vector FASTA",
            key=_FORM_KEYS["vector_fasta"],
            placeholder=r"D:\TLA\vectors\vector.fa  or  /mnt/d/TLA/vectors/vector.fa",
            help="Vector/GOI FASTA file shared across all samples in this run.",
        )
        run_name = st.text_input(
            "Run name",
            key=_FORM_KEYS["run_name"],
            placeholder="e.g. PCSK9_exp1",
            help=(
                "Short name for this run (no slashes). Output will be saved to "
                f"`{OUTPUT_ROOT}/tla-pip/<run_name>/`."
            ),
        )
        if run_name.strip():
            st.caption(
                f"Output \u2192 `{OUTPUT_ROOT}/tla-pip/{run_name.strip()}/`"
            )
        st.subheader("Optional inputs")

        references_dir = st.text_input(
            "References directory (optional)",
            key=_FORM_KEYS["references_dir"],
            placeholder="Leave blank to auto-download into output directory",
            help=(
                "Local host-genome and mixed-ref cache directory. "
                "When blank the pipeline downloads and caches inside the output directory."
            ),
        )
        accession = st.text_input(
            "Host genome NCBI accession",
            key=_FORM_KEYS["accession"],
            placeholder=_DEFAULT_ACCESSION,
            help=(
                f"NCBI accession for the host genome. "
                f"Default: `{_DEFAULT_ACCESSION}` (CHO-K1, *Cricetulus griseus*)."
            ),
        )
        force_download = st.checkbox(
            "Force re-download of host genome",
            value=False,
            help="Re-download even if a cached copy already exists.",
        )

        submitted = st.form_submit_button("Submit")

    if not submitted:
        return None

    return {
        "sample_sheet": sample_sheet,
        "fastq_dir": fastq_dir,
        "vector_fasta": vector_fasta,
        "run_name": run_name,
        "references_dir": references_dir,
        "accession": accession,
        "force_download": force_download,
    }


def parse_command(cmd_str: str) -> dict:
    """Parse a tla-pip command string into a prefill dict for render_form().

    Maps CLI flag names (hyphenated) to form field names.  output_dir is
    always cleared so the user must explicitly choose a new destination.
    """
    from utils import parse_command_flags
    flags = parse_command_flags(cmd_str)
    # Best-effort: extract run name from the tail of any --output-dir value.
    out_raw = flags.get("output-dir", "")
    run_name = str(pathlib.Path(out_raw).name) if out_raw else ""
    return {
        "sample_sheet":   flags.get("sample-sheet", ""),
        "fastq_dir":      flags.get("fastq-dir", ""),
        "vector_fasta":   flags.get("vector-fasta", ""),
        "run_name":       run_name,
        "references_dir": flags.get("references-dir", ""),
        "accession":      flags.get("accession", _DEFAULT_ACCESSION),
    }


def validate_and_build_command(
    raw_values: dict,
) -> tuple[dict, list[str], str]:
    """Validate and normalise submitted values, then build the pipeline command.

    Returns:
        normalized_config : dict   Values saved with the job.
        command           : list   Pipeline command as an argument list.
        working_directory : str    Directory from which the command should run.

    Raises ValueError with a user-facing message on invalid input.
    """
    # --- Required fields --------------------------------------------------
    for field, label in [
        ("sample_sheet", "Sample sheet"),
        ("fastq_dir", "FASTQ directory"),
        ("vector_fasta", "Vector FASTA"),
        ("run_name", "Run name"),
    ]:
        if not raw_values.get(field, "").strip():
            raise ValueError(f"{label} is required.")

    run_name = raw_values["run_name"].strip()
    if any(c in run_name for c in ("/", "\\")):
        raise ValueError("Run name must not contain path separators (/ or \\).")

    # --- Path normalisation -----------------------------------------------
    sample_sheet = normalize_path(raw_values["sample_sheet"])
    fastq_dir = normalize_path(raw_values["fastq_dir"])
    vector_fasta = normalize_path(raw_values["vector_fasta"])
    output_dir = build_output_dir(PIPELINE_ID, run_name)

    references_dir_raw = raw_values.get("references_dir", "").strip()
    references_dir = normalize_path(references_dir_raw) if references_dir_raw else None

    # --- Existence checks -------------------------------------------------
    if not sample_sheet.is_file():
        raise ValueError(f"Sample sheet not found: {sample_sheet}")
    if not fastq_dir.is_dir():
        raise ValueError(f"FASTQ directory not found: {fastq_dir}")
    if not vector_fasta.is_file():
        raise ValueError(f"Vector FASTA not found: {vector_fasta}")
    if references_dir is not None and not references_dir.is_dir():
        raise ValueError(f"References directory not found: {references_dir}")

    # --- Accession --------------------------------------------------------
    accession = raw_values.get("accession", "").strip() or _DEFAULT_ACCESSION

    # --- Entry point check ------------------------------------------------
    if not _PIPELINE_SCRIPT.is_file():
        raise ValueError(
            f"TLA pipeline script not found: {_PIPELINE_SCRIPT}\n"
            "Is the tla-pip directory present?"
        )

    python_exe = _tla_python()

    # --- Build command ----------------------------------------------------
    command: list[str] = [
        python_exe, "-u", str(_PIPELINE_SCRIPT),
        "--sample-sheet", str(sample_sheet),
        "--fastq-dir", str(fastq_dir),
        "--vector-fasta", str(vector_fasta),
        "--output-dir", str(output_dir),
        "--accession", accession,
    ]
    if references_dir is not None:
        command += ["--references-dir", str(references_dir)]
    if raw_values.get("force_download"):
        command.append("--force-download")

    normalized_config = {
        "sample_sheet": str(sample_sheet),
        "fastq_dir": str(fastq_dir),
        "vector_fasta": str(vector_fasta),
        "run_name": run_name,
        "output_dir": str(output_dir),
        "references_dir": str(references_dir) if references_dir else None,
        "accession": accession,
        "force_download": bool(raw_values.get("force_download")),
    }

    return normalized_config, command, str(_PIPELINE_DIR)
