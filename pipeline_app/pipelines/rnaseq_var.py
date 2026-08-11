"""RNAseq Variant Calling pipeline module for the pipeline Streamlit app.

This wraps **RNAseq-var** (BWA-based, general RNA-seq SNP/indel calling).
It is NOT the RNAseq-pip chimeric-transcript pipeline.

Entry point : RNAseq-var/pipeline_run.py
Working dir : <run_dir>  (pipeline reads sample_xlsx/*.xlsx relative to CWD)

Input files supplied by the user
---------------------------------
    FASTQ directory    Flat directory of paired FASTQ files.
    Sample sheet       *.xlsx with fastq (forward/reverse), GOI reference file,
                       and optional Clone ID columns.
    Reference FASTA    Host reference genome FASTA.
    GOI FASTA          Gene-of-interest FASTA.

At submission time the app runs organize_inputs.py to symlink these into the
required subdirectory layout under the run directory:

    <OUTPUT_ROOT>/rnaseq-var/<run_name>/
        Fastq/          <- symlinks to all FASTQ files in --fastq-dir
        reference/      <- symlink to reference FASTA
        goi/            <- symlink to GOI FASTA
        sample_xlsx/    <- symlink to sample sheet
        Result/         <- pipeline writes outputs here

Sample sheet columns (required)
--------------------------------
    fastq (forward)     R1 FASTQ filename (must be in Fastq/)
    fastq (reverse)     R2 FASTQ filename (must be in Fastq/)
    GOI reference file  GOI FASTA filename (must match the file in goi/)
    Clone ID (optional) Replicate-aware SNP grouping when filled for ALL samples.
"""

import pathlib
import subprocess

import streamlit as st

from utils import normalize_path, build_output_dir, OUTPUT_ROOT, collect_mixed_output_files

PIPELINE_ID  = "rnaseq-var"
DISPLAY_NAME = "RNAseq Variant Calling (RNAseq-var)"

_CONDA_ENV         = "RNAseq"
_PIPELINE_DIR      = pathlib.Path(__file__).parents[2] / "RNAseq-var"
_PIPELINE_SCRIPT   = _PIPELINE_DIR / "pipeline_run.py"
_ORGANIZE_SCRIPT   = _PIPELINE_DIR / "organize_inputs.py"

# Widget keys for form pre-fill via st.session_state.
_FORM_KEYS = {
    "fastq_dir":    "rnaseqvar_fastq_dir",
    "sample_sheet": "rnaseqvar_sample_sheet",
    "reference":    "rnaseqvar_reference",
    "goi":          "rnaseqvar_goi",
    "run_name":     "rnaseqvar_run_name",
}


def collect_export_files(output_dir):
    """Return final result files for export, excluding intermediates.

    RNAseq-var writes intermediates (BAMs, trimmed FASTQs, BWA indexes,
    reference FASTA copies) alongside final results in the same Result/
    directory tree.  This delegates to collect_mixed_output_files, which
    walks the full tree and strips known intermediate types by suffix.

    Final files kept:
      snp.txt, *.fastp.json, *.fastp.html, read_counts.txt,
      *.coverage_plot.png, *.vcf (snp / indel / consensus / amino-annotated),
      pindel output files (non-BAM/FASTQ).
    """
    return collect_mixed_output_files(output_dir)


def _rnaseq_python() -> str:
    """Return the Python executable from the RNAseq conda environment.

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
            f"Create it first: conda env create -f RNAseq-var/RNAseq.yml\n"
            f"{result.stderr.strip()}"
        )
    return result.stdout.strip()


def render_form(prefill: dict | None = None) -> dict | None:
    """Render the RNAseq Variant Calling submission form.

    Returns submitted raw values as a dict, or None before submission.
    """
    if prefill:
        for field, key in _FORM_KEYS.items():
            val = prefill.get(field)
            if val is not None:
                st.session_state[key] = val

    with st.form("rnaseqvar_form"):
        st.info(
            "**RNAseq Variant Calling (RNAseq-var)** — BWA-based SNP and indel calling "
            "from RNA-seq data for any gene of interest.  \n"
            "This is **not** the chimeric-transcript pipeline (RNAseq-pip)."
        )

        st.subheader("Input files")

        fastq_dir = st.text_input(
            "FASTQ directory",
            key=_FORM_KEYS["fastq_dir"],
            placeholder=r"Z:\RNAseq-var\run1\Fastq  or  /mnt/z/RNAseq-var/run1/Fastq",
            help="Flat directory containing all paired FASTQ files for this run.",
        )
        sample_sheet = st.text_input(
            "Sample sheet (.xlsx)",
            key=_FORM_KEYS["sample_sheet"],
            placeholder=r"Z:\RNAseq-var\run1\samples.xlsx",
            help=(
                "Excel file with columns: fastq (forward), fastq (reverse), "
                "GOI reference file, Clone ID (optional)."
            ),
        )
        reference = st.text_input(
            "Reference genome FASTA",
            key=_FORM_KEYS["reference"],
            placeholder=r"Z:\refs\GRCh38_latest_genomic.fna",
            help="Host reference genome FASTA used for BWA alignment.",
        )
        goi = st.text_input(
            "GOI FASTA",
            key=_FORM_KEYS["goi"],
            placeholder=r"Z:\refs\BCMA_LTR.fa",
            help="Gene-of-interest FASTA used for targeted variant calling.",
        )

        st.subheader("Run name")
        run_name = st.text_input(
            "Run name",
            key=_FORM_KEYS["run_name"],
            placeholder="e.g. F8_RNAseq_var_run1",
            help=(
                "Short label for this run (no slashes). "
                f"Output will be saved to `{OUTPUT_ROOT}/rnaseq-var/<run_name>/`."
            ),
        )
        if run_name.strip():
            st.caption(
                f"Output -> `{OUTPUT_ROOT}/rnaseq-var/{run_name.strip()}/Result/`"
            )

        submitted = st.form_submit_button("Submit")

    if not submitted:
        return None

    return {
        "fastq_dir":    fastq_dir,
        "sample_sheet": sample_sheet,
        "reference":    reference,
        "goi":          goi,
        "run_name":     run_name,
    }


def parse_command(cmd_str: str) -> dict:
    """Parse a previous pipeline_run.py command string into a prefill dict."""
    # pipeline_run.py takes a single positional argument (the run_dir).
    # The 4 source input paths are not in the command, so leave them blank.
    parts = cmd_str.strip().split()
    run_dir = ""
    for part in reversed(parts):
        if not part.startswith("-"):
            run_dir = part
            break
    return {
        "fastq_dir":    "",
        "sample_sheet": "",
        "reference":    "",
        "goi":          "",
        "run_name":     pathlib.Path(run_dir).name if run_dir else "",
    }


def validate_and_build_command(
    raw_values: dict,
) -> tuple[dict, list[str], str]:
    """Validate inputs, organize the run directory, and build the pipeline command.

    organize_inputs.py is run synchronously at submission time to create the
    required subdirectory structure (symlinks only, no copy).  The tmux job
    then runs pipeline_run.py against the already-organized run directory.

    Returns:
        normalized_config : dict   Values saved with the job history entry.
        command           : list   pipeline_run.py command as an argument list.
        working_directory : str    run_dir (pipeline reads sample_xlsx/*.xlsx
                                   relative to CWD).

    Raises ValueError with a user-facing message on invalid input.
    """
    for field, label in [
        ("fastq_dir",    "FASTQ directory"),
        ("sample_sheet", "Sample sheet"),
        ("reference",    "Reference FASTA"),
        ("goi",          "GOI FASTA"),
        ("run_name",     "Run name"),
    ]:
        if not raw_values.get(field, "").strip():
            raise ValueError(f"{label} is required.")

    run_name = raw_values["run_name"].strip()
    if any(c in run_name for c in ("/", "\\")):
        raise ValueError("Run name must not contain path separators (/ or \\).")

    fastq_dir    = normalize_path(raw_values["fastq_dir"])
    sample_sheet = normalize_path(raw_values["sample_sheet"])
    reference    = normalize_path(raw_values["reference"])
    goi          = normalize_path(raw_values["goi"])
    run_dir      = build_output_dir(PIPELINE_ID, run_name)

    if not fastq_dir.is_dir():
        raise ValueError(f"FASTQ directory not found: {fastq_dir}")
    if not sample_sheet.is_file():
        raise ValueError(f"Sample sheet not found: {sample_sheet}")
    if sample_sheet.suffix.lower() != ".xlsx":
        raise ValueError(f"Sample sheet must be an .xlsx file: {sample_sheet}")
    if not reference.is_file():
        raise ValueError(f"Reference FASTA not found: {reference}")
    if not goi.is_file():
        raise ValueError(f"GOI FASTA not found: {goi}")

    if not _PIPELINE_SCRIPT.is_file():
        raise ValueError(
            f"RNAseq-var pipeline script not found: {_PIPELINE_SCRIPT}\n"
            "Is the RNAseq-var directory present?"
        )
    if not _ORGANIZE_SCRIPT.is_file():
        raise ValueError(
            f"organize_inputs.py not found: {_ORGANIZE_SCRIPT}"
        )

    python_exe = _rnaseq_python()

    # --- Organize inputs synchronously before the tmux job is launched -------
    organize_result = subprocess.run(
        [
            python_exe, str(_ORGANIZE_SCRIPT),
            "--fastq-dir",    str(fastq_dir),
            "--sample-sheet", str(sample_sheet),
            "--reference",    str(reference),
            "--goi",          str(goi),
            "--run-dir",      str(run_dir),
            "--overwrite",
        ],
        capture_output=True,
        text=True,
    )
    if organize_result.returncode != 0:
        raise ValueError(
            "Input organization failed:\n"
            + organize_result.stdout
            + organize_result.stderr
        )

    # --- Pipeline command for tmux -------------------------------------------
    command: list[str] = [
        python_exe, "-u", str(_PIPELINE_SCRIPT),
        str(run_dir),
    ]

    # The pipeline writes to <run_dir>/Result/.
    output_dir = str(run_dir / "Result")

    normalized_config = {
        "fastq_dir":    str(fastq_dir),
        "sample_sheet": str(sample_sheet),
        "reference":    str(reference),
        "goi":          str(goi),
        "run_name":     run_name,
        "output_dir":   output_dir,
    }

    # working_directory = run_dir so glob("sample_xlsx/*.xlsx") in
    # pipeline_run.py resolves correctly.
    return normalized_config, command, str(run_dir)
