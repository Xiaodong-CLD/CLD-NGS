"""Amplicon sequencing pipeline module for the pipeline Streamlit app.

Entry point : amp-pip-dev/run_pipeline.py
Working dir : amp-pip-dev/

Required CLI args:
    --sample-sheet   Sample sheet (.xlsx)
    --fastq-dir      Directory containing paired FASTQ files
    --reference      Amplicon reference FASTA file
    --output-dir     Output directory (created automatically)

The pipeline expects the sample sheet to have these columns:
    fastq (forward)         R1 filename (must match a file in fastq-dir)
    fastq (reverse)         R2 filename (must match a file in fastq-dir)
    GOI reference file      reference FASTA filename (must match reference)
    Clone ID (optional)     replicate grouping key (leave blank to disable grouping)

Output layout after the run:

    <output_dir>/
        Results/
            snp.txt                     run-level SNP report (variants ≥ 2% VarFreq)
            indel.txt                   run-level Pindel structural-variant report
            <sample>/
                *.fastp.html / .json    fastp QC reports
                coverage-plot-output/   per-sample coverage plots
                varscan-output/         SNP, indel, and amino-annotated VCFs
                pindel-output/          pindel-output.Summary.vcf (ad_frac ≥ 1%)
        tmp/
            <sample>/                   all pipeline intermediates (BAMs, mpileup,
                                        trimmed FASTQs, pindel intermediates).
                                        Safe to delete after reviewing Results/.
"""

import pathlib
import subprocess

import streamlit as st

from utils import normalize_path, build_output_dir, OUTPUT_ROOT

PIPELINE_ID  = "amp-pip"
DISPLAY_NAME = "Amplicon Pipeline"

_CONDA_ENV    = "variant-calling"
_PIPELINE_DIR = pathlib.Path(__file__).parents[2] / "amp-pip-dev"
_PIPELINE_SCRIPT = _PIPELINE_DIR / "pipeline_combined.py"

# Widget keys for form pre-fill via st.session_state.
_FORM_KEYS = {
    "sample_sheet": "amp_sample_sheet",
    "fastq_dir":    "amp_fastq_dir",
    "reference":    "amp_reference",
    "run_name":     "amp_run_name",
}


_FASTQ_EXTS = (".fastq", ".fastq.gz", ".fq", ".fq.gz")


def _prepare_input_dir(
    fastq_dir: pathlib.Path,
    sample_sheet: pathlib.Path,
    reference: pathlib.Path,
    output_dir: pathlib.Path,
) -> None:
    """Create the symlink layout that pipeline_combined.py expects.

    output_dir/
        Fastq/          symlinks to every FASTQ in fastq_dir
        sample_xlsx/    symlink to sample_sheet
        reference/      symlink to reference

    Idempotent: re-running updates existing symlinks (ln -sfn semantics).
    Raises ValueError if no FASTQ files are found.
    """
    fastq_target = output_dir / "Fastq"
    sheet_target = output_dir / "sample_xlsx"
    ref_target   = output_dir / "reference"
    for d in (fastq_target, sheet_target, ref_target):
        d.mkdir(parents=True, exist_ok=True)

    # Symlink FASTQ files.
    n_fq = 0
    for src in sorted(fastq_dir.iterdir()):
        if any(src.name.endswith(ext) for ext in _FASTQ_EXTS):
            dest = fastq_target / src.name
            if dest.exists() and not dest.is_symlink():
                raise ValueError(f"Regular file already occupies FASTQ slot: {dest}")
            if dest.is_symlink():
                dest.unlink()
            dest.symlink_to(src.resolve())
            n_fq += 1
    if n_fq == 0:
        raise ValueError(
            f"No FASTQ files (.fastq/.fastq.gz/.fq/.fq.gz) found in: {fastq_dir}"
        )

    # Symlink sample sheet.
    sheet_dest = sheet_target / sample_sheet.name
    if sheet_dest.exists() and not sheet_dest.is_symlink():
        raise ValueError(f"Regular file already occupies sample sheet slot: {sheet_dest}")
    if sheet_dest.is_symlink():
        sheet_dest.unlink()
    sheet_dest.symlink_to(sample_sheet.resolve())

    # Symlink reference FASTA.
    ref_dest = ref_target / reference.name
    if ref_dest.exists() and not ref_dest.is_symlink():
        raise ValueError(f"Regular file already occupies reference slot: {ref_dest}")
    if ref_dest.is_symlink():
        ref_dest.unlink()
    ref_dest.symlink_to(reference.resolve())


def _amp_python() -> str:
    """Return the Python executable from the variant-calling conda environment.

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
            f"Create it first: conda env create -f amp-pip-dev/envs/variant-calling.yml\n"
            f"{result.stderr.strip()}"
        )
    return result.stdout.strip()


def render_form(prefill: dict | None = None) -> dict | None:
    """Render the Amplicon Pipeline submission form.

    Returns submitted raw values as a dict, or None before submission.
    """
    if prefill:
        for field, key in _FORM_KEYS.items():
            val = prefill.get(field)
            if val is not None:
                st.session_state[key] = val

    with st.form("amp_form"):
        st.subheader("Required inputs")

        sample_sheet = st.text_input(
            "Sample sheet (.xlsx)",
            key=_FORM_KEYS["sample_sheet"],
            placeholder=r"Z:\Amplicon\run1\samples.xlsx  or  /mnt/z/Amplicon/run1/samples.xlsx",
            help="Path to the sample sheet. Windows and WSL paths are both accepted.",
        )
        fastq_dir = st.text_input(
            "FASTQ directory",
            key=_FORM_KEYS["fastq_dir"],
            placeholder=r"Z:\Amplicon\run1\Fastq  or  /mnt/z/Amplicon/run1/Fastq",
            help="Directory containing the paired FASTQ files listed in the sample sheet.",
        )
        reference = st.text_input(
            "Reference FASTA",
            key=_FORM_KEYS["reference"],
            placeholder=r"Z:\Amplicon\refs\goi.fa  or  /mnt/z/Amplicon/refs/goi.fa",
            help="Amplicon reference FASTA file (GOI).",
        )
        run_name = st.text_input(
            "Run name",
            key=_FORM_KEYS["run_name"],
            placeholder="e.g. PCSK9_amplicon_run1",
            help=(
                "Short label for this run (no slashes). Output will be saved to "
                f"`{OUTPUT_ROOT}/amp-pip/<run_name>/`."
            ),
        )
        if run_name.strip():
            st.caption(
                f"Output \u2192 `{OUTPUT_ROOT}/amp-pip/{run_name.strip()}/`"
            )

        submitted = st.form_submit_button("Submit")

    if not submitted:
        return None

    return {
        "sample_sheet": sample_sheet,
        "fastq_dir":    fastq_dir,
        "reference":    reference,
        "run_name":     run_name,
    }


def parse_command(cmd_str: str) -> dict:
    """Parse a previous run_pipeline.py command string into a prefill dict."""
    from utils import parse_command_flags
    flags = parse_command_flags(cmd_str)
    out_raw = flags.get("output-dir", "")
    run_name = str(pathlib.Path(out_raw).name) if out_raw else ""
    return {
        "sample_sheet": flags.get("sample-sheet", ""),
        "fastq_dir":    flags.get("fastq-dir", ""),
        "reference":    flags.get("reference", ""),
        "run_name":     run_name,
    }


def validate_and_build_command(
    raw_values: dict,
) -> tuple[dict, list[str], str]:
    """Validate submitted values and build the run_pipeline.py command.

    Returns:
        normalized_config : dict   Values saved with the job history entry.
        command           : list   Command as an argument list.
        working_directory : str    Directory from which the command should run.

    Raises ValueError with a user-facing message on invalid input.
    """
    for field, label in [
        ("sample_sheet", "Sample sheet"),
        ("fastq_dir",    "FASTQ directory"),
        ("reference",    "Reference FASTA"),
        ("run_name",     "Run name"),
    ]:
        if not raw_values.get(field, "").strip():
            raise ValueError(f"{label} is required.")

    run_name = raw_values["run_name"].strip()
    if any(c in run_name for c in ("/", "\\")):
        raise ValueError("Run name must not contain path separators (/ or \\).")

    sample_sheet = normalize_path(raw_values["sample_sheet"])
    fastq_dir    = normalize_path(raw_values["fastq_dir"])
    reference    = normalize_path(raw_values["reference"])
    output_dir   = build_output_dir(PIPELINE_ID, run_name)

    if not sample_sheet.is_file():
        raise ValueError(f"Sample sheet not found: {sample_sheet}")
    if not fastq_dir.is_dir():
        raise ValueError(f"FASTQ directory not found: {fastq_dir}")
    if not reference.is_file():
        raise ValueError(f"Reference FASTA not found: {reference}")

    if not _PIPELINE_SCRIPT.is_file():
        raise ValueError(
            f"Amplicon pipeline script not found: {_PIPELINE_SCRIPT}\n"
            "Is the amp-pip-dev directory present?"
        )

    # Prepare the input directory layout expected by pipeline_combined.py.
    _prepare_input_dir(fastq_dir, sample_sheet, reference, output_dir)

    python_exe = _amp_python()

    command: list[str] = [
        python_exe, "-u", str(_PIPELINE_SCRIPT),
        str(output_dir),
    ]

    normalized_config = {
        "sample_sheet": str(sample_sheet),
        "fastq_dir":    str(fastq_dir),
        "reference":    str(reference),
        "run_name":     run_name,
        "output_dir":   str(output_dir),
    }

    return normalized_config, command, str(_PIPELINE_DIR)
