"""KO pipeline module for the pipeline Streamlit app.

Entry point : ko-pip/pipeline.py
Working dir : ko-pip/

Required CLI args:
    -i / --input_dir    Directory containing paired FASTQ files
    -s / --sample_sheet Path to the sample sheet (.xlsx)
    -r / --reference    Reference FASTA
    --blast_ref         Directory of plasmid FASTAs (*.fa/*.fasta/*.fna) for BLAST

Optional CLI args:
    -o / --output_dir   Output directory (default: current working directory)
"""

import pathlib
import subprocess

import streamlit as st

from utils import normalize_path, build_output_dir, OUTPUT_ROOT

PIPELINE_ID = "ko-pip"
DISPLAY_NAME = "Knockout Pipeline"

_CONDA_ENV = "ko_pip"
_PIPELINE_DIR = pathlib.Path(__file__).parents[2] / "ko-pip"
_PIPELINE_SCRIPT = _PIPELINE_DIR / "pipeline.py"

# Widget keys for form pre-fill via st.session_state.
_FORM_KEYS = {
    "fastq_dir":    "ko_fastq_dir",
    "sample_sheet": "ko_sample_sheet",
    "reference":    "ko_reference",
    "blast_ref":    "ko_blast_ref",
    "run_name":     "ko_run_name",
}


def _ko_python() -> str:
    """Return the Python executable from the ko_pip conda environment.

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
            f"Create it first: conda env create -f ko-pip/envs/ko_pip.yml\n"
            f"{result.stderr.strip()}"
        )
    return result.stdout.strip()


def render_form(prefill: dict | None = None) -> dict | None:
    """Render the KO submission form.

    If *prefill* is provided (e.g. from parsing a previous command), the form
    widgets are pre-populated with those values so the user only needs to
    change the output directory before submitting.

    Returns submitted raw values as a dict, or None before submission.
    """
    # Pre-fill widget session state before the form renders.  Setting
    # st.session_state[key] before st.text_input(key=key) is the official
    # Streamlit way to programmatically populate a widget.
    if prefill:
        for field, key in _FORM_KEYS.items():
            val = prefill.get(field)
            if val is not None:
                st.session_state[key] = val

    with st.form("ko_form"):
        st.subheader("Required inputs")

        fastq_dir = st.text_input(
            "FASTQ directory",
            key=_FORM_KEYS["fastq_dir"],
            placeholder=r"Z:\KO\run1\Fastq  or  /mnt/z/KO/run1/Fastq",
            help="Directory containing the paired-end FASTQ files for this run. Windows and WSL paths are both accepted.",
        )
        sample_sheet = st.text_input(
            "Sample sheet (.xlsx)",
            key=_FORM_KEYS["sample_sheet"],
            placeholder=r"Z:\KO\run1\samples.xlsx  or  /mnt/z/KO/run1/samples.xlsx",
            help="Full path to the sample sheet Excel file.",
        )
        reference = st.text_input(
            "Reference FASTA",
            key=_FORM_KEYS["reference"],
            placeholder=r"Z:\KO\refs\reference.fa  or  /mnt/z/KO/refs/reference.fa",
            help="Reference genome FASTA for alignment.",
        )
        blast_ref = st.text_input(
            "BLAST reference directory",
            key=_FORM_KEYS["blast_ref"],
            placeholder=r"Z:\KO\blast_ref  or  /mnt/z/KO/blast_ref",
            help=(
                "Directory of plasmid FASTA files (*.fa / *.fasta / *.fna) "
                "used to classify breakpoint soft-clip clusters."
            ),
        )

        st.subheader("Run name")

        run_name = st.text_input(
            "Run name",
            key=_FORM_KEYS["run_name"],
            placeholder="e.g. F8KO_pilot",
            help=(
                "Short name for this run (no slashes). Output will be saved to "
                f"`{OUTPUT_ROOT}/ko-pip/<run_name>/`."
            ),
        )
        if run_name.strip():
            st.caption(
                f"Output → `{OUTPUT_ROOT}/ko-pip/{run_name.strip()}/`"
            )

        submitted = st.form_submit_button("Submit")

    if not submitted:
        return None

    return {
        "fastq_dir": fastq_dir,
        "sample_sheet": sample_sheet,
        "reference": reference,
        "blast_ref": blast_ref,
        "run_name": run_name,
    }


def parse_command(cmd_str: str) -> dict:
    """Parse a ko-pip command string into a prefill dict for render_form().

    Maps CLI flag names to form field names.  output_dir is always cleared
    so the user must explicitly choose a new destination.
    """
    from utils import parse_command_flags
    flags = parse_command_flags(cmd_str)
    # Best-effort: extract run name from the tail of any --output_dir value.
    out_raw = flags.get("output_dir", "")
    run_name = str(pathlib.Path(out_raw).name) if out_raw else ""
    return {
        "fastq_dir":    flags.get("input_dir", ""),
        "sample_sheet": flags.get("sample_sheet", ""),
        "reference":    flags.get("reference", ""),
        "blast_ref":    flags.get("blast_ref", ""),
        "run_name":     run_name,
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
        ("fastq_dir",     "FASTQ directory"),
        ("sample_sheet",  "Sample sheet"),
        ("reference",     "Reference FASTA"),
        ("blast_ref",     "BLAST reference directory"),
        ("run_name",      "Run name"),
    ]:
        if not raw_values.get(field, "").strip():
            raise ValueError(f"{label} is required.")

    run_name = raw_values["run_name"].strip()
    if any(c in run_name for c in ("/", "\\")):
        raise ValueError("Run name must not contain path separators (/ or \\).")

    # --- Path normalisation -----------------------------------------------
    fastq_dir    = normalize_path(raw_values["fastq_dir"])
    sample_sheet = normalize_path(raw_values["sample_sheet"])
    reference    = normalize_path(raw_values["reference"])
    blast_ref    = normalize_path(raw_values["blast_ref"])
    output_dir   = build_output_dir(PIPELINE_ID, run_name)

    # --- Existence checks -------------------------------------------------
    if not fastq_dir.is_dir():
        raise ValueError(f"FASTQ directory not found: {fastq_dir}")

    if not sample_sheet.is_file():
        raise ValueError(f"Sample sheet not found: {sample_sheet}")
    if sample_sheet.suffix.lower() != ".xlsx":
        raise ValueError(f"Sample sheet must be an .xlsx file: {sample_sheet}")

    if not reference.is_file():
        raise ValueError(f"Reference FASTA not found: {reference}")

    if not blast_ref.is_dir():
        raise ValueError(f"BLAST reference directory not found: {blast_ref}")

    # --- Entry point check ------------------------------------------------
    if not _PIPELINE_SCRIPT.is_file():
        raise ValueError(
            f"KO pipeline script not found: {_PIPELINE_SCRIPT}\n"
            "Is the ko-pip directory present?"
        )

    python_exe = _ko_python()

    # --- Build command ----------------------------------------------------
    command: list[str] = [
        python_exe, "-u", str(_PIPELINE_SCRIPT),
        "--input_dir",    str(fastq_dir),
        "--sample_sheet", str(sample_sheet),
        "--reference",    str(reference),
        "--blast_ref",    str(blast_ref),
        "--output_dir",   str(output_dir),
    ]

    normalized_config = {
        "fastq_dir":    str(fastq_dir),
        "sample_sheet": str(sample_sheet),
        "reference":    str(reference),
        "blast_ref":    str(blast_ref),
        "run_name":     run_name,
        "output_dir":   str(output_dir),
    }

    return normalized_config, command, str(_PIPELINE_DIR)
