"""Shared submission utilities for the pipeline Streamlit app."""

import json
import logging
import os
import pathlib
import re
import shlex
import shutil
import subprocess
import tempfile
import time
import traceback
import uuid
import zipfile
from datetime import datetime, timezone

# ---------------------------------------------------------------------------
# Windows Z: drive → Linux path mapping  (PERSONAL — adjust for your setup)
#
# On this server the Z: drive is mounted at /mnt/z.  The Windows path root
# Z:\ maps directly to /mnt/z — there is NO extra sub-path offset.
# Concretely:
#
#   Windows Z:\  ==  Linux /mnt/z/
#
# So "Z:\CLD\MiSeq\MiSeqOutput\run42\Fastq"
#   →  "/mnt/z/CLD/MiSeq/MiSeqOutput/run42/Fastq"
#
# How to find YOUR own values if this ever changes:
#   1. On the Linux server run:  mount | grep cifs   (or: df -h)
#      The output shows the actual mount point for the share.
#   2. Set _Z_MOUNT_ROOT to that mount point (e.g. /mnt/z).
#   3. Set _WIN_Z_ROOT to the sub-directory WITHIN the mount that corresponds
#      to the Z:\ root on Windows.  If Z:\ maps directly to /mnt/z with no
#      offset (the current setup), use pathlib.Path("") (empty string).
#
# Tip — how to copy an absolute Windows path without any ambiguity:
#   • In File Explorer: click the address bar to reveal the full path
#     (e.g. Z:\CLD\MiSeq\run42\Fastq), then Ctrl+C.
#   • Or: Shift + Right-click on the folder → "Copy as path" (gives a
#     quoted absolute path you can paste directly into this app).
#   • Both forms (forward- or back-slashes, with or without quotes) are
#     handled by normalize_path() below.
# ---------------------------------------------------------------------------
_Z_MOUNT_ROOT = pathlib.Path("/mnt/z")
_WIN_Z_ROOT = pathlib.Path("")   # <-- empty = Z:\ maps directly to /mnt/z/
                                  #     set to e.g. "CLD/MiSeq" if your mount
                                  #     has an extra sub-path offset

# ---------------------------------------------------------------------------
# Disk-space and stale-cleanup thresholds  (adjust as needed)
# ---------------------------------------------------------------------------
STALE_DAYS = 30              # tmp/ dirs older than this (by mtime) are eligible for removal
Z_DRIVE_LOW_SPACE_GB = 10.0  # warn when Z: drive free space falls below this
LOCAL_LOW_SPACE_GB = 5.0    # warn when local pipeline-output filesystem falls below this


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def _strip_quotes(s: str) -> str:
    """Strip surrounding whitespace and quote characters from a path string.

    Handles:
      - Straight single/double quotes  ' "  (e.g. from Windows "Copy as path")
      - Unicode curly/smart quotes  \u2018\u2019 \u201c\u201d
    Strips repeatedly until no more matching pairs are found.
    """
    _QUOTE_PAIRS = (
        ('"', '"'),
        ("'", "'"),
        ("\u2018", "\u2019"),   # ' '
        ("\u201c", "\u201d"),   # " "
    )
    s = s.strip()
    changed = True
    while changed:
        changed = False
        for open_q, close_q in _QUOTE_PAIRS:
            if len(s) >= 2 and s[0] == open_q and s[-1] == close_q:
                s = s[1:-1].strip()
                changed = True
                break
    return s


def _is_windows_path(s: str) -> bool:
    """Return True if *s* looks like a Windows drive path (e.g. C:\\... or C:/...)."""
    return len(s) >= 3 and s[1] == ":" and s[2] in ("/", "\\")


def _z_mount_available(drive_letter: str) -> bool:
    """Return True if /mnt/z is mounted and accessible.

    Only the Z: drive is checked since that is the configured Windows share mount.
    Other drive letters fall back to wslpath.
    """
    if drive_letter.upper() != "Z":
        return False
    return _Z_MOUNT_ROOT.is_dir()


# ---------------------------------------------------------------------------
# Z: drive health check & remount
# ---------------------------------------------------------------------------

def check_z_drive() -> dict:
    """Return a dict describing the Z: drive mount state.

    Keys:
      mounted   bool        – True when /mnt/z exists and is stat-able
      writable  bool        – True when a test file can be created (only if mounted)
      free_gb   float|None  – available space in GiB (None if stat failed)
      total_gb  float|None  – total capacity in GiB (None if stat failed)
      error     str         – non-empty when an unexpected error occurred
    """
    result = {"mounted": False, "writable": False,
              "free_gb": None, "total_gb": None, "error": ""}
    try:
        if not _Z_MOUNT_ROOT.is_dir():
            return result
        result["mounted"] = True
        # Disk-space probe.
        try:
            usage = shutil.disk_usage(_Z_MOUNT_ROOT)
            result["free_gb"] = usage.free / (1024 ** 3)
            result["total_gb"] = usage.total / (1024 ** 3)
        except OSError:
            pass
        # Quick writable check: create a zero-byte probe file and remove it.
        probe = _Z_MOUNT_ROOT / ".pipeline_app_probe"
        try:
            probe.touch()
            probe.unlink()
            result["writable"] = True
        except OSError:
            result["writable"] = False
    except Exception as exc:
        result["error"] = str(exc)
    return result


def check_local_disk() -> dict:
    """Return disk-space info for the filesystem that hosts OUTPUT_ROOT.

    Keys:
      free_gb   float|None  – available space in GiB (None if stat failed)
      total_gb  float|None  – total capacity in GiB
      error     str         – non-empty when stat failed
    """
    result: dict = {"free_gb": None, "total_gb": None, "error": ""}
    try:
        OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
        usage = shutil.disk_usage(OUTPUT_ROOT)
        result["free_gb"] = usage.free / (1024 ** 3)
        result["total_gb"] = usage.total / (1024 ** 3)
    except OSError as exc:
        result["error"] = str(exc)
    return result



def _windows_to_z_path(s: str) -> str:
    """Convert a Z: Windows path to its Linux equivalent using the configured offsets.

    E.g. "Z:\\CLD\\MiSeq\\MiSeqOutput\\run42\\Fastq"
         ->  "/mnt/z/CLD/MiSeq/MiSeqOutput/run42/Fastq"

    The mapping is controlled by _Z_MOUNT_ROOT and _WIN_Z_ROOT at the top of
    this file — adjust those two constants to match your server's mount layout.
    """
    rest = s[3:].replace("\\", "/")   # strip "Z:\" or "Z:/"
    return str(_Z_MOUNT_ROOT / _WIN_Z_ROOT / rest)


def normalize_path(raw: str) -> pathlib.Path:
    """Convert a user-supplied path (Windows or WSL/Linux) to a pathlib.Path.

    Conversion order for Windows paths:
      1. If /z/<drive>/ is mounted on this server, map to /z/<drive>/...
      2. Otherwise fall back to ``wslpath -u`` (standard WSL /mnt/<drive>/ layout).

    Linux/WSL paths are left unchanged.

    Raises ValueError if conversion fails.
    """
    s = _strip_quotes(raw)
    if not s:
        raise ValueError("Path is empty.")

    if _is_windows_path(s):
        drive = s[0]
        if _z_mount_available(drive):
            s = _windows_to_z_path(s)
        else:
            result = subprocess.run(
                ["wslpath", "-u", s],
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                raise ValueError(
                    f"Could not convert Windows path {s!r}: {result.stderr.strip()}\n"
                    f"/mnt/z is not mounted and wslpath also failed."
                )
            s = result.stdout.strip()

    return pathlib.Path(s)


# ---------------------------------------------------------------------------
# Output root & run history
# ---------------------------------------------------------------------------

# Resolve relative to this file so every repo clone writes to its own
# pipeline_outputs/ directory (ngs/pipeline_outputs/) regardless of CWD
# or who is running the app.
OUTPUT_ROOT = pathlib.Path(__file__).parents[1] / "pipeline_outputs"

# Run-history (job metadata + logs) lives inside the output root.
_RUNS_DIR = OUTPUT_ROOT / "run_history"


def build_output_dir(pipeline_id: str, run_name: str) -> pathlib.Path:
    """Return the resolved output path for a run: OUTPUT_ROOT/pipeline_id/run_name.

    Creates the pipeline sub-directory (not the run folder itself — the
    pipeline does that) so the parent always exists when the job starts.
    """
    out = OUTPUT_ROOT / pipeline_id / run_name
    out.parent.mkdir(parents=True, exist_ok=True)
    return out


def _make_job_id(pipeline_id: str) -> str:
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
    short_uid = uuid.uuid4().hex[:6]
    return f"{ts}_{pipeline_id}_{short_uid}"


def _run_dir(job_id: str) -> pathlib.Path:
    parts = job_id.split("_")
    # parts: [timestamp, pipeline_id, uid]  (pipeline_id may itself contain _)
    # group: date portion (first 8 chars of timestamp)
    date_part = parts[0][:8]  # YYYYMMDD
    pipeline_id = parts[1]
    return _RUNS_DIR / pipeline_id / date_part / job_id


# ---------------------------------------------------------------------------
# launch.sh generation
# ---------------------------------------------------------------------------

def _write_launch_sh(
    run_dir: pathlib.Path,
    working_dir: str,
    command: list[str],
    status_file: pathlib.Path,
    log_file: pathlib.Path,
    exit_code_file: pathlib.Path,
) -> pathlib.Path:
    """Write a shell script that runs the pipeline and records status/exit code."""
    launch_path = run_dir / "launch.sh"

    # Safely quote all command arguments.
    quoted_cmd = shlex.join(command)
    status_path = shlex.quote(str(status_file))
    log_path = shlex.quote(str(log_file))
    exit_code_path = shlex.quote(str(exit_code_file))

    # If command[0] is an absolute path (e.g. a conda env's Python), prepend
    # its bin/ directory to PATH so sibling tools (fastp, bwa, samtools, etc.)
    # from the same environment are found without a full conda activate.
    path_export = ""
    exe = command[0] if command else ""
    if exe.startswith("/"):
        env_bin = shlex.quote(str(pathlib.Path(exe).parent))
        path_export = f"export PATH={env_bin}:$PATH\n"

    script = f"""\
#!/bin/bash
set -euo pipefail
cd {shlex.quote(working_dir)}
{path_export}printf 'running' > {status_path}
{quoted_cmd} >> {log_path} 2>&1
_exit=$?
printf '%d' $_exit > {exit_code_path}
if [ $_exit -eq 0 ]; then
    printf 'completed' > {status_path}
else
    printf 'failed' > {status_path}
fi
exit $_exit
"""
    launch_path.write_text(script)
    launch_path.chmod(0o755)
    return launch_path


# ---------------------------------------------------------------------------
# Job submission
# ---------------------------------------------------------------------------

def submit_job(
    pipeline_id: str,
    display_name: str,
    normalized_config: dict,
    command: list[str],
    working_directory: str,
    output_dir: str,
) -> dict:
    """Create run directory, write metadata files, and start a detached tmux session.

    Returns a job metadata dict (same content as job.json).
    Callers are responsible for checking check_output_dir_conflict() first.
    """
    job_id = _make_job_id(pipeline_id)
    run_dir = _run_dir(job_id)
    run_dir.mkdir(parents=True, exist_ok=True)

    session_name = f"pipeline_{pipeline_id}_{job_id[-6:]}"
    submitted_at = datetime.now(timezone.utc).isoformat()
    human_cmd = shlex.join(command)

    log_file = run_dir / "pipeline.log"
    status_file = run_dir / "status.json"
    exit_code_file = run_dir / "exit_code.txt"
    command_file = run_dir / "command.txt"

    # Write readable command.
    command_file.write_text(human_cmd + "\n")

    # Write initial status.
    status_file.write_text("submitted")

    # Create empty log so tail works immediately.
    log_file.touch()

    # Write launch script.
    launch_path = _write_launch_sh(
        run_dir, working_directory, command,
        status_file, log_file, exit_code_file,
    )

    # Metadata.
    job_meta = {
        "job_id": job_id,
        "pipeline_id": pipeline_id,
        "display_name": display_name,
        "submitted_at": submitted_at,
        "config": normalized_config,
        "output_dir": output_dir,
        "command": command,
        "human_command": human_cmd,
        "working_directory": working_directory,
        "session_name": session_name,
        "run_dir": str(run_dir),
        "log_file": str(log_file),
    }
    (run_dir / "job.json").write_text(
        json.dumps(job_meta, indent=2) + "\n"
    )

    # Start detached tmux session.
    subprocess.run(
        [
            "tmux", "new-session", "-d",
            "-s", session_name,
            str(launch_path),
        ],
        check=True,
    )

    return job_meta


def check_output_dir_conflict(pipeline_id: str, output_dir: str) -> dict | None:
    """Return a conflict info dict if *output_dir* is associated with a previous job.

    Returns a dict with keys: job_id, status, exit_code (int or None).
    Returns None when no previous job used this output directory.

    Active statuses (submitted, running) mean the job is still in progress.
    Terminal statuses (completed, failed, interrupted) are informational warnings.
    """
    pipeline_runs = _RUNS_DIR / pipeline_id
    if not pipeline_runs.exists():
        return None
    for job_json in sorted(pipeline_runs.rglob("job.json"),
                           key=lambda p: p.stat().st_mtime, reverse=True):
        try:
            meta = json.loads(job_json.read_text())
        except (json.JSONDecodeError, OSError):
            continue
        if meta.get("output_dir") == output_dir:
            status = read_job_status(meta)
            exit_code = None
            ec_file = pathlib.Path(meta["run_dir"]) / "exit_code.txt"
            if ec_file.exists():
                try:
                    exit_code = int(ec_file.read_text().strip())
                except (ValueError, OSError):
                    pass
            return {"job_id": meta["job_id"], "status": status, "exit_code": exit_code}
    return None


# ---------------------------------------------------------------------------
# Job status
# ---------------------------------------------------------------------------

def _tmux_session_exists(session_name: str) -> bool:
    result = subprocess.run(
        ["tmux", "has-session", "-t", session_name],
        capture_output=True,
    )
    return result.returncode == 0


def read_job_status(job_meta: dict) -> str:
    """Return one of: submitted, running, completed, failed, interrupted, unknown."""
    run_dir = pathlib.Path(job_meta["run_dir"])
    status_file = run_dir / "status.json"
    exit_code_file = run_dir / "exit_code.txt"

    if not status_file.exists():
        return "unknown"

    try:
        raw_status = status_file.read_text().strip()
    except OSError:
        return "unknown"

    if raw_status in ("completed", "failed"):
        return raw_status

    if raw_status == "running":
        if exit_code_file.exists():
            # Exit code written but status not updated — treat as completed/failed.
            try:
                code = int(exit_code_file.read_text().strip())
                return "completed" if code == 0 else "failed"
            except (ValueError, OSError):
                pass
        # No exit code yet — check tmux.
        if _tmux_session_exists(job_meta["session_name"]):
            return "running"
        # Session gone without exit code → interrupted.
        return "interrupted"

    if raw_status == "submitted":
        return "submitted"

    return "unknown"


def kill_job(job_meta: dict) -> str:
    """Kill the tmux session for a job and mark it interrupted.

    Returns a short human-readable result string:
      - 'killed'       session was running and has been killed
      - 'not_running'  session was already gone
      - 'error'        tmux kill-session failed
    """
    session = job_meta["session_name"]
    run_dir = pathlib.Path(job_meta["run_dir"])
    status_file = run_dir / "status.json"

    if not _tmux_session_exists(session):
        return "not_running"

    result = subprocess.run(
        ["tmux", "kill-session", "-t", session],
        capture_output=True,
    )
    if result.returncode != 0:
        return "error"

    # Mark status as interrupted so the UI reflects the kill immediately.
    try:
        status_file.write_text("interrupted")
    except OSError:
        pass
    return "killed"


def delete_job(job_meta: dict) -> None:
    """Remove the run directory (metadata + log only — pipeline output is untouched)."""
    run_dir = pathlib.Path(job_meta["run_dir"])
    if run_dir.exists():
        shutil.rmtree(run_dir)
        # Prune empty parent dirs (date bucket, pipeline bucket).
        for parent in [run_dir.parent, run_dir.parent.parent]:
            try:
                parent.rmdir()   # no-op unless the directory is now empty
            except OSError:
                break


# ---------------------------------------------------------------------------
# Stale temporary-file cleanup
# ---------------------------------------------------------------------------

def cleanup_stale_tmp(stale_days: int = STALE_DAYS) -> list[str]:
    """Remove ``tmp/`` directories from pipeline output dirs of old terminal jobs.

    A ``tmp/`` directory at ``<output_dir>/tmp/`` is removed when ALL of the
    following are true:

    1. Its job appears in run_history with a terminal status
       (completed, failed, or interrupted) — never active jobs.
    2. The job's ``submitted_at`` timestamp is older than *stale_days* days.
       (Filesystem mtime is NOT used — ``tmp/`` mtime is set at pipeline
       start and is unreliable for measuring how long a job has been done.)
    3. The ``tmp/`` directory is under OUTPUT_ROOT (safety guard).

    Nothing outside OUTPUT_ROOT is touched.  Removal errors are logged but
    do not raise.  Returns a list of paths successfully removed (as strings).
    """
    log = _export_logger()
    removed: list[str] = []

    if not _RUNS_DIR.exists():
        return removed

    now_ts = datetime.now(timezone.utc).timestamp()
    cutoff_ts = now_ts - stale_days * 86400
    output_root_resolved = OUTPUT_ROOT.resolve()

    for jf in _RUNS_DIR.rglob("job.json"):
        try:
            meta = json.loads(jf.read_text())
        except (json.JSONDecodeError, OSError):
            continue

        # Only clean terminal jobs.
        status = read_job_status(meta)
        if status not in ("completed", "failed", "interrupted"):
            continue

        # Age gate: use submitted_at, not filesystem mtime.
        try:
            submitted_ts = datetime.fromisoformat(meta["submitted_at"]).timestamp()
        except (KeyError, ValueError):
            continue
        if submitted_ts >= cutoff_ts:
            continue

        output_dir = meta.get("output_dir", "")
        if not output_dir:
            continue
        tmp_dir = pathlib.Path(output_dir) / "tmp"
        if not tmp_dir.is_dir():
            continue

        # Safety: tmp_dir must be under OUTPUT_ROOT.
        try:
            tmp_dir.resolve().relative_to(output_root_resolved)
        except ValueError:
            continue

        age_days = (now_ts - submitted_ts) / 86400
        try:
            shutil.rmtree(tmp_dir)
            removed.append(str(tmp_dir))
            log.info(
                "Removed stale tmp dir (job %.0f days old): %s",
                age_days, tmp_dir,
            )
        except OSError as exc:
            log.warning("Could not remove stale tmp dir %s: %s", tmp_dir, exc)

    return removed


# ---------------------------------------------------------------------------
# Results copy
# ---------------------------------------------------------------------------

def copy_results(output_dir: str, dest_raw: str) -> pathlib.Path:
    """Copy *output_dir* to *dest_raw* (Linux or Windows Z: path), excluding any
    ``tmp`` sub-directory.  Returns the resolved destination path.

    Uses shutil.copytree with ``dirs_exist_ok=True`` so it is safe to call
    repeatedly (e.g. to update a partially-copied destination).
    """
    src = pathlib.Path(output_dir)
    if not src.exists():
        raise ValueError(f"Output directory not found: {src}")
    dst = normalize_path(dest_raw)
    dst.mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        src, dst,
        ignore=shutil.ignore_patterns("tmp"),
        dirs_exist_ok=True,
    )
    return dst


# ---------------------------------------------------------------------------
# Export final results to the Windows desktop  (WSL only)
#
# Layout on the desktop:
#   <Desktop>/CLD_NGS_pipeline_outputs_temp/<pipeline_id>/<output_name>.zip
#
# Only final, user-relevant outputs are packaged.  Large intermediate files
# (BAM/SAM/CRAM, FASTQ, pileups, reference indexes, caller working dirs, the
# tmp/ directory, ...) are excluded via an explicit allow/deny list — the ZIP
# is built directly from the selected files, never "zip everything then delete".
# ---------------------------------------------------------------------------

EXPORT_FOLDER_NAME = "CLD_NGS_pipeline_outputs_temp"

# Top-level directories inside a run's output_dir that hold final results.
# Matched case-insensitively.
_INCLUDE_TOP_DIRS = {
    "final_report", "final_reports", "final_results",
    "results", "report", "reports", "summary",
    "logs", "figures", "plots", "qc",
}

# Directory names that are always excluded, anywhere in the tree
# (temporary files, reference builds/indexes, caller working dirs, caches).
_EXCLUDE_DIR_NAMES = {
    "tmp", "temp", "reference", "__pycache__",
    ".cache", ".snakemake", "work", "workdir",
}

# Top-level file suffixes considered final/user-relevant (README, tables,
# configs, run parameters, small figures, interpretation docs).  Matched
# case-insensitively.
_INCLUDE_TOP_FILE_SUFFIXES = {
    ".tsv", ".csv", ".txt", ".json", ".yml", ".yaml",
    ".md", ".pdf", ".html", ".xlsx", ".png", ".jpg", ".jpeg",
    ".log", ".fa", ".fasta",
}

# File suffixes that are always excluded, anywhere in the tree — large or
# intermediate alignment/sequence/index/cache artifacts not needed for
# result interpretation.  Checked with a case-insensitive endswith().
_EXCLUDE_FILE_SUFFIXES = (
    ".bam", ".sam", ".cram", ".bai", ".crai", ".csi", ".tbi", ".idx",
    ".mpileup", ".pileup",
    ".fastq", ".fq", ".fastq.gz", ".fq.gz",
    ".sa", ".pac", ".bwt", ".ann", ".amb", ".dict", ".fai",
    ".bt2", ".bt2l",
    ".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto", ".njs",
)

# Characters that are invalid in Windows filenames.
_WIN_INVALID_CHARS = set('<>:"/\\|?*')
_WIN_RESERVED_NAMES = {
    "CON", "PRN", "AUX", "NUL",
    *(f"COM{i}" for i in range(1, 10)),
    *(f"LPT{i}" for i in range(1, 10)),
}


class ExportError(Exception):
    """User-facing export failure. Technical details are written to the app log."""


class ExportZipExists(ExportError):
    """Raised when the target ZIP already exists and overwrite was not requested."""

    def __init__(self, zip_path: pathlib.Path, windows_path: str):
        self.zip_path = zip_path
        self.windows_path = windows_path
        super().__init__(f"A ZIP named '{zip_path.name}' already exists.")


def _export_logger() -> logging.Logger:
    """Return a file-backed logger for export diagnostics (created lazily)."""
    lg = logging.getLogger("pipeline_app.export")
    if not lg.handlers:
        lg.setLevel(logging.INFO)
        log_path = OUTPUT_ROOT / "run_history" / "app.log"
        try:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            handler = logging.FileHandler(log_path)
            handler.setFormatter(
                logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
            )
            lg.addHandler(handler)
        except OSError:
            lg.addHandler(logging.NullHandler())
    return lg


def log_app_exception(context: str) -> None:
    """Record the current exception traceback in the app log (concise UI aside)."""
    _export_logger().error("%s:\n%s", context, traceback.format_exc())


def sanitize_output_name(name: str) -> str:
    """Return *name* with only Windows-invalid characters replaced.

    Replaces reserved characters (``< > : " / \\ | ? *``) and control
    characters with ``_``.  Trailing dots/spaces (invalid on Windows) are
    trimmed, reserved device names are prefixed, and a trailing ``.zip`` the
    user may have typed is dropped.  Falls back to ``results`` when empty.
    The name is otherwise left unchanged.
    """
    name = (name or "").strip()
    if name.lower().endswith(".zip"):
        name = name[:-4]
    cleaned = "".join(
        "_" if (c in _WIN_INVALID_CHARS or ord(c) < 32) else c
        for c in name
    )
    cleaned = cleaned.rstrip(" .")
    if not cleaned:
        return "results"
    if cleaned.upper() in _WIN_RESERVED_NAMES:
        cleaned = f"_{cleaned}"
    return cleaned


def _run_windows_tool(args: list[str]) -> str | None:
    """Run a Windows interop tool (reg.exe / cmd.exe), trying System32 fallback.

    Returns stdout (str) on success, or None if the tool is unavailable or fails.
    """
    for candidate in (args[0], f"/mnt/c/Windows/System32/{args[0]}"):
        try:
            result = subprocess.run(
                [candidate, *args[1:]],
                capture_output=True, text=True, timeout=15,
            )
        except (FileNotFoundError, subprocess.SubprocessError):
            continue
        if result.returncode == 0:
            return result.stdout.replace("\r", "")
    return None


def _expand_win_env(value: str) -> str:
    """Expand %VAR% tokens in a Windows path via cmd.exe (best effort)."""
    if "%" not in value:
        return value
    out = _run_windows_tool(["cmd.exe", "/c", "echo", value])
    return out.strip() if out else value


def _query_desktop_registry() -> str:
    """Return the resolved Windows Desktop path from the registry, or "".

    Reads ``Shell Folders\\Desktop`` which reflects OneDrive redirection when
    the Desktop has been moved.
    """
    out = _run_windows_tool([
        "reg.exe", "query",
        r"HKCU\Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders",
        "/v", "Desktop",
    ])
    if not out:
        return ""
    for line in out.splitlines():
        line = line.strip()
        if line.lower().startswith("desktop"):
            parts = re.split(r"\s+REG_(?:EXPAND_)?SZ\s+", line, maxsplit=1)
            if len(parts) == 2 and parts[1].strip():
                return _expand_win_env(parts[1].strip())
    return ""


def _wslpath(mode: str, path: str) -> str | None:
    """Convert a path with ``wslpath`` (mode '-u' to WSL, '-w' to Windows)."""
    try:
        result = subprocess.run(
            ["wslpath", mode, path],
            capture_output=True, text=True, timeout=15,
        )
    except (FileNotFoundError, subprocess.SubprocessError):
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip()


def detect_windows_desktop() -> dict:
    """Detect the Windows Desktop path and its WSL equivalent (dynamic, no hard-coding).

    Returns a dict with keys:
      windows_path  str   – e.g. ``C:\\Users\\Robolab\\Desktop`` ("" if unknown)
      wsl_path      str   – e.g. ``/mnt/c/Users/Robolab/Desktop`` ("" if unknown)
      writable      bool  – True when a probe file can be created there
      onedrive      bool  – True when the Desktop is redirected into OneDrive
      error         str   – non-empty on failure
    """
    info = {"windows_path": "", "wsl_path": "", "writable": False,
            "onedrive": False, "error": ""}
    try:
        win_path = _query_desktop_registry()
        if not win_path:
            profile = _run_windows_tool(["cmd.exe", "/c", "echo", "%USERPROFILE%"])
            profile = (profile or "").strip()
            if profile and not profile.startswith("%"):
                win_path = profile.rstrip("\\") + "\\Desktop"
        if not win_path:
            info["error"] = "Could not determine the Windows Desktop path."
            return info
        info["windows_path"] = win_path
        info["onedrive"] = "onedrive" in win_path.lower()

        wsl_path = _wslpath("-u", win_path)
        if not wsl_path:
            info["error"] = f"Could not convert Desktop path to WSL: {win_path}"
            return info
        info["wsl_path"] = wsl_path

        desktop = pathlib.Path(wsl_path)
        if not desktop.is_dir():
            info["error"] = f"Desktop path does not exist: {wsl_path}"
            return info
        probe = desktop / ".pipeline_app_probe"
        try:
            probe.touch()
            probe.unlink()
            info["writable"] = True
        except OSError:
            info["writable"] = False
        return info
    except Exception as exc:  # pragma: no cover - defensive
        _export_logger().error("detect_windows_desktop failed:\n%s", traceback.format_exc())
        info["error"] = str(exc)
        return info


# Additional FASTA-like suffixes excluded only from the mixed-output collector.
# These are reference copies written by pipelines into their result directory
# (not user-facing outputs) that are NOT in _EXCLUDE_FILE_SUFFIXES because
# the clean-output pipelines (amp, ko, tla) never write them to Results/.
_EXCLUDE_MIXED_REF_SUFFIXES = (".fasta", ".fa", ".fna", ".fsa")


def _is_excluded_file(name: str) -> bool:
    lower = name.lower()
    return any(lower.endswith(suf) for suf in _EXCLUDE_FILE_SUFFIXES)


def collect_final_files(output_dir: str | os.PathLike) -> list[tuple[pathlib.Path, pathlib.Path]]:
    """Return the final, user-relevant files under *output_dir* to include in a ZIP.

    Each item is ``(absolute_path, relative_path)`` where *relative_path* is
    relative to *output_dir*.  Selection rules:

      • Top-level result directories in _INCLUDE_TOP_DIRS are walked; within
        them every file is kept unless it lives in an excluded sub-directory
        (_EXCLUDE_DIR_NAMES) or has an excluded suffix (_EXCLUDE_FILE_SUFFIXES).
      • Top-level files whose suffix is in _INCLUDE_TOP_FILE_SUFFIXES (or named
        README*) are kept.
      • Everything else (tmp/, reference/, caller working dirs, BAM/FASTQ/...)
        is skipped.
    """
    root = pathlib.Path(output_dir)
    selected: list[tuple[pathlib.Path, pathlib.Path]] = []
    if not root.is_dir():
        return selected

    for entry in sorted(root.iterdir()):
        if entry.is_dir():
            if entry.name.lower() in _EXCLUDE_DIR_NAMES:
                continue
            if entry.name.lower() not in _INCLUDE_TOP_DIRS:
                continue
            for path in sorted(entry.rglob("*")):
                if not path.is_file():
                    continue
                rel = path.relative_to(root)
                if any(part.lower() in _EXCLUDE_DIR_NAMES for part in rel.parts[:-1]):
                    continue
                if _is_excluded_file(path.name):
                    continue
                selected.append((path, rel))
        elif entry.is_file():
            suffix = entry.suffix.lower()
            is_readme = entry.stem.lower() == "readme"
            if (suffix in _INCLUDE_TOP_FILE_SUFFIXES or is_readme) \
                    and not _is_excluded_file(entry.name):
                selected.append((entry, entry.relative_to(root)))
    return selected


def collect_mixed_output_files(
    output_dir: str | os.PathLike,
) -> list[tuple[pathlib.Path, pathlib.Path]]:
    """Return final files from a pipeline where intermediates and results are mixed.

    Unlike :func:`collect_final_files`, which uses an allowlist of known result
    directory names, this function walks the *entire* output_dir tree and
    **excludes** known intermediate types.  Use this for pipelines like
    RNAseq-var where VCFs, QC reports, BAMs, and FASTQ intermediates all live
    in the same directory tree.

    Excluded:
      * Directories in _EXCLUDE_DIR_NAMES (tmp/, reference/, etc.)
      * File suffixes in _EXCLUDE_FILE_SUFFIXES (BAM, FASTQ, BWA indexes, ...)
      * Reference FASTA copies (_EXCLUDE_MIXED_REF_SUFFIXES: .fasta .fa .fna .fsa)
    """
    root = pathlib.Path(output_dir)
    selected: list[tuple[pathlib.Path, pathlib.Path]] = []
    if not root.is_dir():
        return selected

    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(root)
        # Skip if any ancestor directory is on the exclude list.
        if any(part.lower() in _EXCLUDE_DIR_NAMES for part in rel.parts[:-1]):
            continue
        # Skip standard intermediate suffixes (BAM, FASTQ, index files, ...).
        if _is_excluded_file(path.name):
            continue
        # Skip reference FASTA copies that pipelines write alongside results.
        lower = path.name.lower()
        if any(lower.endswith(suf) for suf in _EXCLUDE_MIXED_REF_SUFFIXES):
            continue
        selected.append((path, rel))

    return selected


def export_results_zip(
    output_dir: str,
    pipeline_id: str,
    output_name: str,
    desktop_wsl_path: str,
    overwrite: bool = False,
    file_collector=None,
) -> dict:
    """Package the final results of a run into a ZIP on the Windows desktop.

    Creates ``<desktop>/CLD_NGS_pipeline_outputs_temp/<pipeline_id>/`` and
    writes ``<output_name>.zip`` (name sanitized for Windows) containing only
    the files selected by *file_collector* (or :func:`collect_final_files` when
    *file_collector* is None), each stored under a top-level ``<output_name>/``
    folder inside the archive.

    *file_collector* — optional callable ``(output_dir: Path) -> list[(Path, Path)]``.
    Pass a pipeline-specific collector (e.g. :func:`collect_mixed_output_files`)
    for pipelines that mix intermediates and final files in the same directory.
    When None, :func:`collect_final_files` is used (the default for all pipelines
    whose outputs are already separated into Results/ vs tmp/).

    Raises :class:`ExportZipExists` if the ZIP already exists and *overwrite*
    is False.  Raises :class:`ExportError` for missing results, an unwritable
    destination, insufficient disk space, or a ZIP/copy failure.  Full
    tracebacks are recorded in the app log; the exception messages are concise.

    Returns a dict: ``{zip_path, wsl_path, windows_path, size_bytes, file_count}``.
    """
    log = _export_logger()
    safe_name = sanitize_output_name(output_name)

    src = pathlib.Path(output_dir)
    if not src.is_dir():
        raise ExportError(f"Pipeline output directory not found: {output_dir}")

    collector = file_collector if file_collector is not None else collect_final_files
    files = collector(src)
    if not files:
        raise ExportError(
            "No final result files were found to export for this run."
        )

    export_dir = pathlib.Path(desktop_wsl_path) / EXPORT_FOLDER_NAME / pipeline_id
    try:
        export_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        log.error("Could not create export directory %s:\n%s", export_dir, traceback.format_exc())
        raise ExportError(
            f"Could not create the export folder on the desktop: {export_dir}. "
            "Check that the desktop is available and writable."
        ) from exc

    zip_path = export_dir / f"{safe_name}.zip"
    windows_path = _wslpath("-w", str(zip_path)) or str(zip_path)

    if zip_path.exists() and not overwrite:
        raise ExportZipExists(zip_path, windows_path)

    # Disk-space guard: require free space >= total (uncompressed) size.
    total_size = 0
    for abs_path, _ in files:
        try:
            total_size += abs_path.stat().st_size
        except OSError:
            pass
    try:
        free = shutil.disk_usage(export_dir).free
    except OSError:
        free = None
    if free is not None and free < total_size:
        raise ExportError(
            "Not enough free disk space on the desktop to create this ZIP "
            f"(need ~{total_size // (1024 * 1024)} MB, "
            f"free {free // (1024 * 1024)} MB)."
        )

    # Build to a temporary file in the same directory, then atomically replace
    # so a failed run never leaves a partial ZIP in place.
    tmp_fd, tmp_name = tempfile.mkstemp(suffix=".zip.part", dir=str(export_dir))
    os.close(tmp_fd)
    tmp_path = pathlib.Path(tmp_name)
    try:
        with zipfile.ZipFile(tmp_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for abs_path, rel in files:
                zf.write(abs_path, arcname=str(pathlib.Path(safe_name) / rel))
        os.replace(tmp_path, zip_path)
    except OSError as exc:
        try:
            tmp_path.unlink()
        except OSError:
            pass
        log.error("ZIP creation failed for %s:\n%s", zip_path, traceback.format_exc())
        raise ExportError(
            "Failed to create the ZIP file. See the app log for details."
        ) from exc

    size_bytes = zip_path.stat().st_size
    log.info("Exported %d files to %s (%d bytes)", len(files), zip_path, size_bytes)
    return {
        "zip_path": str(zip_path),
        "wsl_path": str(zip_path),
        "windows_path": windows_path,
        "size_bytes": size_bytes,
        "file_count": len(files),
    }


# ---------------------------------------------------------------------------
# Log tail
# ---------------------------------------------------------------------------

_LOG_TAIL_LINES = 150


def tail_log(job_meta: dict, n: int = _LOG_TAIL_LINES) -> str:
    """Return approximately the last *n* lines of the pipeline log."""
    log_path = pathlib.Path(job_meta["log_file"])
    if not log_path.exists():
        return "(log file not found)"

    result = subprocess.run(
        ["tail", "-n", str(n), str(log_path)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return f"(could not read log: {result.stderr.strip()})"
    return result.stdout


# ---------------------------------------------------------------------------
# Recent jobs listing
# ---------------------------------------------------------------------------

def list_recent_jobs(pipeline_id: str | None = None, limit: int = 20) -> list[dict]:
    """Return the most recent *limit* job metadata dicts, newest first."""
    search_root = _RUNS_DIR / pipeline_id if pipeline_id else _RUNS_DIR
    if not search_root.exists():
        return []

    job_files = sorted(
        search_root.rglob("job.json"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    jobs = []
    for jf in job_files[:limit]:
        try:
            jobs.append(json.loads(jf.read_text()))
        except (json.JSONDecodeError, OSError):
            continue
    return jobs


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

def format_submitted_at(iso_str: str) -> str:
    """Return a compact local-time string from an ISO UTC timestamp.

    Example: "2026-07-29  14:57:14.382 BST"
    """
    try:
        dt_utc = datetime.fromisoformat(iso_str)
        dt_local = dt_utc.astimezone()   # converts to system local timezone
        tz_name = dt_local.strftime("%Z")
        return dt_local.strftime(f"%Y-%m-%d  %H:%M:%S.") + f"{dt_local.microsecond // 1000:03d} {tz_name}"
    except (ValueError, TypeError):
        return iso_str


def format_human_command(command: list[str]) -> str:
    """Return a multiline shell representation of the pipeline command.

    Skips the Python executable and '-u' flag (not useful to display).
    Pairs up --flag value arguments across separate lines for readability.

    Example:
        python pipeline.py \\
          --sample-sheet /path/to/samples.xlsx \\
          --fastq-dir /path/to/Fastq \\
          --output-dir /path/to/output
    """
    # Drop leading python/-u so display starts from the script name.
    args = list(command)
    while args and args[0] in ("-u",):
        args.pop(0)
    if args and not args[0].startswith("-"):
        # First non-flag arg is the script path — shorten to basename.
        script = pathlib.Path(args[0]).name
        args[0] = script

    if not args:
        return shlex.join(command)

    lines = [args[0]]
    i = 1
    while i < len(args):
        a = args[i]
        if a.startswith("-") and i + 1 < len(args) and not args[i + 1].startswith("-"):
            lines.append(f"  {a} {shlex.quote(args[i + 1])}")
            i += 2
        else:
            lines.append(f"  {shlex.quote(a)}")
            i += 1

    return " \\\n".join(lines)


# ---------------------------------------------------------------------------
# Command parsing (for "paste previous command" pre-fill)
# ---------------------------------------------------------------------------

def parse_command_flags(cmd_str: str) -> dict[str, str]:
    """Parse --flag value pairs from a shell command string.

    Uses shlex to handle quoted values correctly.  Only long flags (--flag)
    are extracted; single-dash flags like -u are ignored.  A flag immediately
    followed by another --flag (or end of tokens) is treated as a boolean True.

    Returns a dict keyed by flag name (without leading --).
    """
    try:
        tokens = shlex.split(cmd_str)
    except ValueError:
        return {}

    result: dict[str, str | bool] = {}
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok.startswith("--"):
            flag = tok[2:]
            # Next token is the value unless it's another --flag or end.
            if i + 1 < len(tokens) and not tokens[i + 1].startswith("--"):
                result[flag] = tokens[i + 1]
                i += 2
            else:
                result[flag] = True
                i += 1
        else:
            i += 1
    return result
