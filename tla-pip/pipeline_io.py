"""Shared output-integrity helpers for pipeline stages.

- ``is_nonempty``       -- file exists and has size > 0.
- ``gzip_ok``           -- file is a complete, non-truncated gzip/BGZF stream.
- ``bam_ok``            -- file is a complete BAM (valid header + BGZF EOF block).
- ``count_fastq_reads`` -- number of reads in a FASTQ (gzipped or plain).
- ``commit``            -- validate a ``.tmp`` file and atomically rename it
                           into its final path.
- ``write_tsv``         -- write rows + header to a TSV atomically.
- ``get_logger``        -- process-wide timestamped logger (single handler).

Copied from the sibling ko-pip pipeline so the TLA pipeline reuses the same
no-clobber / atomic-write conventions. This module is the shared deposit for
helper functions; it may grow to hold more shared utilities over time.
"""
import os
import gzip
import logging
import subprocess
import sys

_LOGGER_NAME = "tla_pip"


def get_logger(name=_LOGGER_NAME):
    """Return a process-wide configured logger.

    All modules log through this so messages share one timestamped format and
    severity level, with the literal ``[Takeda Pipeline Message]`` tag kept for
    continuity with ko-pip. Output goes to ``stdout`` (not ``stderr``) so the
    tmux wrapper's single redirect captures everything in run order.

    Idempotent: repeated calls reuse the same handler rather than stacking
    duplicate handlers (which would print every line multiple times).
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(stream=sys.stdout)
        handler.setFormatter(logging.Formatter(
            "[%(asctime)s] [Takeda Pipeline Message] [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.propagate = False
    return logger


def is_nonempty(path):
    """True only if ``path`` exists and is non-empty (size > 0)."""
    return os.path.exists(path) and os.path.getsize(path) > 0


def gzip_ok(path):
    """True if ``path`` is a complete, non-truncated gzip/BGZF stream.

    Uses ``gzip -t``, which decompresses the whole stream and fails if the
    end-of-stream marker is missing -- i.e. it detects truncation, not just a
    0-byte file. Works on plain ``.gz`` and on BGZF (``.vcf.gz``) files.
    """
    if not is_nonempty(path):
        return False
    try:
        subprocess.check_call(
            ["gzip", "-t", path],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        return False


def bam_ok(path):
    """True if ``path`` is a complete BAM (valid header + BGZF EOF block).

    ``samtools quickcheck`` is cheap (it seeks; it does not read the whole
    file) and specifically detects a missing BGZF EOF block, which is the
    signature of a truncated BAM.
    """
    if not is_nonempty(path):
        return False
    try:
        subprocess.check_call(
            ["samtools", "quickcheck", path],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        return False


def count_fastq_reads(path):
    """Return the number of reads in a FASTQ file (gzipped or plain).

    Counts newlines in binary chunks (fast) and divides by 4. Raises
    ``ValueError`` if the line count is not a multiple of 4 (truncated record),
    and -- for gzipped input -- the underlying ``gzip`` reader raises if the
    stream itself is truncated. So a successful return means the file is both
    structurally complete and record-aligned.
    """
    opener = gzip.open if path.endswith(".gz") else open
    n_lines = 0
    with opener(path, "rb") as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            n_lines += chunk.count(b"\n")
    if n_lines % 4 != 0:
        raise ValueError(
            f"FASTQ line count {n_lines} is not a multiple of 4 (truncated?): {path}")
    return n_lines // 4


def commit(tmp_path, final_path, validator=is_nonempty):
    """Validate ``tmp_path`` and atomically move it to ``final_path``.

    Returns True on success. On validation failure, removes ``tmp_path`` and
    returns False, leaving ``final_path`` untouched -- so a partial output never
    reaches the canonical name and the next run re-runs the step.
    """
    if not validator(tmp_path):
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        return False
    os.replace(tmp_path, final_path)
    return True


def write_tsv(rows, header, final_path):
    """Write ``rows`` (an iterable of sequences) + ``header`` to a TSV atomically.

    Each row is tab-joined; fields are stringified. The file is written to a
    ``.tmp`` path and renamed into place only after the full write completes,
    so an interrupted write never leaves a partial table at the canonical name.
    """
    os.makedirs(os.path.dirname(os.path.abspath(final_path)), exist_ok=True)
    tmp_path = final_path + ".tmp"
    with open(tmp_path, "w") as fh:
        fh.write("\t".join(str(c) for c in header) + "\n")
        for row in rows:
            fh.write("\t".join(str(c) for c in row) + "\n")
    os.replace(tmp_path, final_path)
    return final_path


if __name__ == "__main__":
    pass
