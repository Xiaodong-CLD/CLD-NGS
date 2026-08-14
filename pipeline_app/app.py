"""Pipeline Launcher – Streamlit app.

Run with:
    cd pipeline_app
    streamlit run app.py
"""

import pathlib

import streamlit as st

import utils as common
from pipelines import tla, ko, amp, rnaseq_var, unmap

# ---------------------------------------------------------------------------
# Pipeline registry – add one entry per new pipeline module.
# ---------------------------------------------------------------------------
PIPELINES = {
    tla.PIPELINE_ID:        tla,
    ko.PIPELINE_ID:         ko,
    amp.PIPELINE_ID:        amp,
    rnaseq_var.PIPELINE_ID: rnaseq_var,
    unmap.PIPELINE_ID:      unmap,
}

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(page_title="CLD NGS Pipeline Launcher", layout="wide")
st.title("CLD NGS Pipeline Launcher")

# ---------------------------------------------------------------------------
# Startup: remove stale tmp/ directories (>30 days, non-active jobs only).
# ---------------------------------------------------------------------------
_cleaned = common.cleanup_stale_tmp()

# Make the kill-job button red via a targeted CSS override.
st.markdown(
    """
    <style>
    button[data-testid="baseButton-primary"][key="kill_most_recent"],
    div[data-testid="column"] button[data-testid="baseButton-primary"]:has(> div > p:contains("Kill")) {
        background-color: #c0392b !important;
        border-color: #c0392b !important;
        color: white !important;
    }
    div[data-testid="stHorizontalBlock"] button[kind="primary"] {
        background-color: #c0392b;
        border-color: #c0392b;
        color: white;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ---------------------------------------------------------------------------
# System checks
# ---------------------------------------------------------------------------

# Z: drive mount status.
_zdrive = common.check_z_drive()
if _zdrive["error"]:
    st.error(f"**Z: drive** — error checking mount: {_zdrive['error']}")
elif _zdrive["mounted"] and _zdrive["writable"]:
    st.success("**Z: drive** — ✅ Mounted and writable")
elif _zdrive["mounted"]:
    st.warning("**Z: drive** — ⚠️ Mounted but **read-only**. Copying results to Z: may fail.")
else:
    st.error("**Z: drive** — ❌ Not mounted. Z: paths will not work.")
    st.warning("Run the following command in a terminal, then click Refresh:")
    st.code("sudo mount -t drvfs Z: /mnt/z", language="bash")
    if st.button("Refresh", key="refresh_after_mount"):
        st.rerun()

# Local Linux disk space (pipeline output storage).
_ldisk = common.check_local_disk()
if _ldisk["error"]:
    st.warning(f"**Local storage** — could not check disk space: {_ldisk['error']}")
elif _ldisk["free_gb"] is not None:
    _l_free = _ldisk["free_gb"]
    _l_total = _ldisk["total_gb"]
    _l_label = (
        f"{_l_free:.1f} / {_l_total:.0f} GB free"
        if _l_total is not None
        else f"{_l_free:.1f} GB free"
    )
    if _l_free < common.LOCAL_LOW_SPACE_GB:
        st.warning(
            f"**Local storage** — ⚠️ Low disk space: {_l_label} on `{common.OUTPUT_ROOT}`"
            f" (threshold: {common.LOCAL_LOW_SPACE_GB:.0f} GB)"
        )
    else:
        st.success(f"**Local storage** — ✅ {_l_label} on `{common.OUTPUT_ROOT}`")

# Report any stale tmp/ directories that were cleaned on startup.
if _cleaned:
    _n = len(_cleaned)
    st.info(
        f"🧹 Removed {_n} stale tmp/ director{'y' if _n == 1 else 'ies'} "
        f"(>{common.STALE_DAYS} days old, non-active jobs only)."
    )

st.divider()

# ---------------------------------------------------------------------------
# Pipeline selector
# ---------------------------------------------------------------------------
pipeline_options = list(PIPELINES.keys())
display_names = {pid: mod.DISPLAY_NAME for pid, mod in PIPELINES.items()}

selected_id = st.selectbox(
    "Select a pipeline",
    options=pipeline_options,
    format_func=lambda pid: display_names[pid],
    key="pipeline_selector",
)

pipeline_mod = PIPELINES[selected_id]

st.divider()


# ---------------------------------------------------------------------------
# Conflict confirmation dialog
# ---------------------------------------------------------------------------
@st.dialog("Output directory already used")
def _confirm_conflict_dialog(conflict: dict, submit_fn):
    """Modal shown when an active job is already using this output directory."""
    status = conflict["status"]
    job_id = conflict["job_id"]
    st.warning(
        f"Job **{job_id}** is currently **{status}** using this output directory. "
        "Submitting again will start a parallel run into the same directory."
    )
    st.markdown("Do you want to submit anyway?")
    col1, col2 = st.columns(2)
    if col1.button("Submit anyway", type="primary"):
        submit_fn()
        st.rerun()
    if col2.button("Cancel"):
        st.rerun()


@st.dialog("Kill job")
def _confirm_kill_dialog(job_meta: dict):
    """Confirmation popup before killing a job."""
    st.warning(
        f"Kill job **{job_meta['job_id']}**?  "
        f"This will terminate tmux session `{job_meta['session_name']}` "
        f"regardless of its current status."
    )
    col1, col2 = st.columns(2)
    if col1.button("Yes, kill it", type="primary"):
        outcome = common.kill_job(job_meta)
        if outcome == "killed":
            st.success("Job killed and marked as interrupted.")
        elif outcome == "not_running":
            st.info("Session was not running (already finished or never started).")
        else:
            st.error("tmux kill-session failed — you may need to kill it manually.")
        st.rerun()
    if col2.button("Cancel"):
        st.rerun()


@st.dialog("Delete run record")
def _confirm_delete_dialog(job_meta: dict):
    """Confirmation popup before deleting a job record from the run history."""
    st.warning(
        f"Delete the record for job **{job_meta['job_id']}**?  \n"
        "This removes the run metadata and log from the app history only.  \n"
        "**The pipeline output directory is not affected.**"
    )
    col1, col2 = st.columns(2)
    if col1.button("Yes, delete", type="primary"):
        common.delete_job(job_meta)
        st.rerun()
    if col2.button("Cancel"):
        st.rerun()


# ---------------------------------------------------------------------------
# Submission helper (called directly or from dialog)
# ---------------------------------------------------------------------------
def _do_submit(pipeline_mod, normalized_config, command, working_dir):
    try:
        common.submit_job(
            pipeline_id=pipeline_mod.PIPELINE_ID,
            display_name=pipeline_mod.DISPLAY_NAME,
            normalized_config=normalized_config,
            command=command,
            working_directory=working_dir,
            output_dir=normalized_config["output_dir"],
        )
        # Signal the paste-command text area to clear on the next render.
        # (Cannot assign to a widget key after it has been instantiated.)
        st.session_state["_reset_paste_cmd"] = True
        st.rerun()
    except Exception as exc:
        st.session_state["_submit_error"] = str(exc)


# ---------------------------------------------------------------------------
# Desktop export (final results ZIP)
# ---------------------------------------------------------------------------
def _do_export(job: dict, output_name: str, overwrite: bool = False) -> None:
    """Run the desktop export and stash the outcome in session_state.

    Never raises; concise messages are stored for display and full tracebacks
    are written to the app log by the utils layer.
    """
    for key in ("_export_result", "_export_error", "_export_conflict"):
        st.session_state.pop(key, None)

    desktop = common.detect_windows_desktop()
    if desktop["error"] or not desktop["wsl_path"]:
        st.session_state["_export_error"] = {
            "job_id": job["job_id"],
            "message": "Could not locate the Windows desktop. "
                       + (desktop["error"] or "unknown reason"),
        }
        return
    if not desktop["writable"]:
        st.session_state["_export_error"] = {
            "job_id": job["job_id"],
            "message": f"The desktop path is not writable: {desktop['windows_path']}",
        }
        return

    try:
        # Let the pipeline module supply its own file collector when it has
        # mixed intermediates/results.  Falls back to the default
        # collect_final_files when the module doesn't define one.
        _pipeline_mod = PIPELINES.get(job.get("pipeline_id", ""))
        _file_collector = getattr(_pipeline_mod, "collect_export_files", None)
        info = common.export_results_zip(
            output_dir=job["output_dir"],
            pipeline_id=job["pipeline_id"],
            output_name=output_name,
            desktop_wsl_path=desktop["wsl_path"],
            overwrite=overwrite,
            file_collector=_file_collector,
        )
        info["job_id"] = job["job_id"]
        st.session_state["_export_result"] = info
    except common.ExportZipExists as exc:
        st.session_state["_export_conflict"] = {
            "job_id": job["job_id"],
            "output_name": output_name,
            "windows_path": exc.windows_path,
        }
    except common.ExportError as exc:
        st.session_state["_export_error"] = {"job_id": job["job_id"], "message": str(exc)}
    except Exception:
        common.log_app_exception("Unexpected error during desktop export")
        st.session_state["_export_error"] = {
            "job_id": job["job_id"],
            "message": "Unexpected error during export. "
                       "Technical details were written to the app log.",
        }


@st.dialog("ZIP already exists")
def _confirm_export_overwrite(job: dict, conflict: dict):
    """Prompt to replace an existing ZIP or cancel and rename."""
    st.warning(
        f"A ZIP for output name **{conflict['output_name']}** already exists:\n\n"
        f"`{conflict['windows_path']}`\n\n"
        "Replace it, or cancel and choose a different output name."
    )
    col1, col2 = st.columns(2)
    if col1.button("Replace existing", type="primary"):
        _do_export(job, conflict["output_name"], overwrite=True)
        st.rerun()
    if col2.button("Cancel"):
        st.session_state.pop("_export_conflict", None)
        st.rerun()


def _render_desktop_export(job: dict) -> None:
    """Render the 'Export final results to desktop' action for a completed job."""
    st.markdown("**Export final results to desktop**")
    st.caption(
        "Packages only the final, user-relevant outputs (reports, tables, "
        "figures, summaries, logs) into a ZIP on the Windows desktop under "
        f"`{common.EXPORT_FOLDER_NAME}\\{job['pipeline_id']}\\`. Large "
        "intermediate files (BAM/FASTQ/pileups/indexes, tmp/) are excluded."
    )

    # Show the outcome of a previous export action for this job.
    result = st.session_state.get("_export_result")
    if result and result.get("job_id") == job["job_id"]:
        size_mb = result["size_bytes"] / (1024 * 1024)
        st.success("Export complete.")
        st.markdown(
            f"- **ZIP filename:** `{pathlib.Path(result['zip_path']).name}`\n"
            f"- **Windows destination:** `{result['windows_path']}`\n"
            f"- **Size:** {size_mb:.2f} MB ({result['file_count']} files)"
        )
        st.warning(
            "This desktop folder is **temporary**. Please move the ZIP to the "
            "shared drive and delete the local copy afterward."
        )
    error = st.session_state.get("_export_error")
    if error and error.get("job_id") == job["job_id"]:
        st.error(error["message"])
    conflict = st.session_state.get("_export_conflict")
    if conflict and conflict.get("job_id") == job["job_id"]:
        _confirm_export_overwrite(job, conflict)

    desktop = common.detect_windows_desktop()
    if desktop["error"] or not desktop["wsl_path"]:
        st.error("Windows desktop not found: " + (desktop["error"] or "unknown reason"))
        return
    if desktop["onedrive"]:
        st.info(f"Desktop is redirected to OneDrive: `{desktop['windows_path']}`")
    if not desktop["writable"]:
        st.error(f"Desktop path is not writable: `{desktop['windows_path']}`")
        return

    default_name = pathlib.Path(job["output_dir"]).name
    name_key = f"export_name_{job['job_id']}"
    out_name = st.text_input(
        "Output name (ZIP filename)",
        value=default_name,
        key=name_key,
        help="Only characters that are invalid in Windows filenames are changed.",
    )
    safe = common.sanitize_output_name(out_name)
    dest = (
        f"{desktop['windows_path']}\\{common.EXPORT_FOLDER_NAME}"
        f"\\{job['pipeline_id']}\\{safe}.zip"
    )
    st.caption(f"Destination: `{dest}`")

    if st.button("Export final results to desktop", key=f"export_btn_{job['job_id']}"):
        with st.spinner("Packaging final results…"):
            _do_export(job, out_name, overwrite=False)
        st.rerun()


# ---------------------------------------------------------------------------
# Submission form (rendered by the pipeline module)
# ---------------------------------------------------------------------------

# "Paste previous command" expander — fills the form from a copied command.
with st.expander("Re-use a previous command (paste to pre-fill form)", expanded=False):
    st.caption(
        "Paste a previous pipeline command below and click **Parse & fill**."
        "All fields except **Run name** will be filled automatically. "
        "Edit as needed, then click Submit."
    )
    # Clear the text area if a submission just completed (flag set by _do_submit).
    if st.session_state.pop("_reset_paste_cmd", False):
        st.session_state["paste_cmd_area"] = ""
    cmd_text = st.text_area(
        "Previous command",
        key="paste_cmd_area",
        height=160,
        placeholder="python -u pipeline.py --input_dir ... --sample_sheet ...",
        label_visibility="collapsed",
    )
    if st.button("Parse & fill form", key="parse_cmd_btn"):
        parse_fn = getattr(pipeline_mod, "parse_command", None)
        if parse_fn is None:
            st.warning("This pipeline does not support command parsing yet.")
        elif not cmd_text.strip():
            st.warning("Paste a command first.")
        else:
            try:
                st.session_state["_prefill"] = parse_fn(cmd_text)
                st.rerun()
            except Exception as exc:
                st.error(f"Could not parse command: {exc}")

_prefill = st.session_state.pop("_prefill", None)
raw_values = pipeline_mod.render_form(prefill=_prefill)

if "_submit_error" in st.session_state:
    st.error(st.session_state.pop("_submit_error"))

if raw_values is not None:
    try:
        normalized_config, command, working_dir = (
            pipeline_mod.validate_and_build_command(raw_values)
        )
    except ValueError as exc:
        st.error(str(exc))
        normalized_config = None

    if normalized_config is not None:
        conflict = common.check_output_dir_conflict(
            pipeline_mod.PIPELINE_ID,
            normalized_config["output_dir"],
        )
        if conflict is not None:
            if conflict["status"] in ("submitted", "running"):
                # Active job — warn before submitting in parallel.
                _confirm_conflict_dialog(
                    conflict,
                    lambda: _do_submit(pipeline_mod, normalized_config, command, working_dir),
                )
            else:
                # Terminal state — silently reuse (resume into same folder).
                _do_submit(pipeline_mod, normalized_config, command, working_dir)
        else:
            _do_submit(pipeline_mod, normalized_config, command, working_dir)


# ---------------------------------------------------------------------------
# Job status panel  (reads from disk — persists across page refreshes)
# ---------------------------------------------------------------------------
st.divider()

_STATUS_COLOR = {
    "submitted": "🟡",
    "running": "🔵",
    "completed": "🟢",
    "failed": "🔴",
    "interrupted": "🟠",
    "unknown": "⚪",
}

recent_jobs = common.list_recent_jobs(pipeline_id=selected_id, limit=10)

if not recent_jobs:
    st.caption("No previous jobs found.")
else:
    # Most recent job — always expanded with Refresh + live log tail.
    most_recent = recent_jobs[0]
    _status = common.read_job_status(most_recent)
    _icon = _STATUS_COLOR.get(_status, "⚪")

    st.subheader(f"Most recent job — `{most_recent['job_id']}`")

    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"**Submitted:** {common.format_submitted_at(most_recent['submitted_at'])}")
        st.markdown(f"**Status:** {_icon} `{_status}`")
        st.markdown(f"**Output dir:** `{most_recent['output_dir']}`")
        st.markdown(f"**tmux session:** `{most_recent['session_name']}`")
    with col2:
        st.markdown(f"**Run dir:** `{most_recent['run_dir']}`")

    st.markdown("**Command:**")
    st.code(
        most_recent.get("human_command") or common.format_human_command(most_recent["command"]),
        language="bash",
    )

    btn_col1, btn_col2, btn_col3 = st.columns([1, 3, 3])
    with btn_col1:
        if st.button("Refresh"):
            st.rerun()
    with btn_col2:
        if st.button("Kill job", type="primary", key="kill_most_recent",
                     help="Terminate the tmux session for this job"):
            _confirm_kill_dialog(most_recent)
    with btn_col3:
        if _status not in ("submitted", "running"):
            if st.button("Delete record", key="delete_most_recent",
                         help="Remove this run record from history (output dir is kept)"):
                _confirm_delete_dialog(most_recent)

    st.markdown("**Recent log output:**")
    _log = common.tail_log(most_recent)
    if _log.strip():
        st.code(_log, language="text")
    else:
        st.caption("(log is empty — pipeline may not have started yet)")

    if _status == "completed":
        st.divider()
        _render_desktop_export(most_recent)
        st.divider()
        st.markdown("**Save results to another location**")
        st.caption(
            "Copies the output directory (excluding `tmp/`) to a destination you choose.  "
            "Accepts Linux paths or Windows Z: paths (e.g. `Z:\\CLD\\NGS-results\\my-run`).  "
            "Re-running copies new or updated files over an existing destination."
        )
        _dest_key = f"dest_{most_recent['job_id']}"
        _dest = st.text_input(
            "Destination path",
            placeholder="Z:\\CLD\\NGS-results\\my-run  or  /mnt/z/CLD/NGS-results/my-run",
            key=_dest_key,
        )
        if st.button("Copy results (excluding tmp)", key=f"copy_{most_recent['job_id']}"):
            if not _dest.strip():
                st.error("Enter a destination path first.")
            else:
                with st.spinner("Copying — this may take a while for large result sets…"):
                    try:
                        resolved = common.copy_results(most_recent["output_dir"], _dest)
                        st.success(f"Results copied to `{resolved}`")
                    except Exception as _exc:
                        st.error(str(_exc))

    # Older jobs — collapsed.
    if len(recent_jobs) > 1:
        st.divider()
        st.subheader(f"Previous {display_names[selected_id]} jobs")
        for _job in recent_jobs[1:]:
            _s = common.read_job_status(_job)
            _lbl = (
                f"{_STATUS_COLOR.get(_s, '⚪')} `{_job['job_id']}`  —  "
                f"{common.format_submitted_at(_job['submitted_at'])}  —  "
                f"**{_s}**"
            )
            with st.expander(_lbl, expanded=False):
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown(f"**Output dir:** `{_job['output_dir']}`")
                    st.markdown(f"**tmux session:** `{_job['session_name']}`")
                with col2:
                    st.markdown(f"**Run dir:** `{_job['run_dir']}`")
                st.markdown("**Command:**")
                st.code(
                    _job.get("human_command") or common.format_human_command(_job["command"]),
                    language="bash",
                )
                st.markdown("**Recent log output:**")
                _log = common.tail_log(_job)
                if _log.strip():
                    st.code(_log, language="text")
                else:
                    st.caption("(log is empty)")
                if _s == "completed":
                    _dest2 = st.text_input(
                        "Copy results to",
                        placeholder="Z:\\CLD\\NGS-results\\my-run  or  /mnt/z/CLD/...",
                        key=f"dest_{_job['job_id']}",
                    )
                    if st.button("Copy results (excluding tmp)", key=f"copy_{_job['job_id']}"):
                        if not _dest2.strip():
                            st.error("Enter a destination path first.")
                        else:
                            with st.spinner("Copying…"):
                                try:
                                    resolved = common.copy_results(_job["output_dir"], _dest2)
                                    st.success(f"Results copied to `{resolved}`")
                                except Exception as _exc:
                                    st.error(str(_exc))
                if _s not in ("submitted", "running"):
                    if st.button("Delete record", key=f"delete_{_job['job_id']}",
                                 help="Remove from history (output dir is kept)"):
                        _confirm_delete_dialog(_job)
