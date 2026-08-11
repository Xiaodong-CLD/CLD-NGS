#!/usr/bin/env python3

"""Build the mixed (host + vector) reference for the TLA pipeline.

The host genome is *not* user-provided: it is downloaded once from NCBI
(default: Chinese hamster CriGri-PICRH-1.0 ``GCF_003668045.3``) via the
``datasets`` CLI, cached on disk, and reused across every run. The user only
supplies the vector / GOI FASTA via ``--vector-fasta`` on the pipeline CLI.

Two on-disk caches, both keyed so work is never repeated:

* ``<host_dir>/<accession>/`` -- the downloaded host genome
  (``host.fa`` + ``.fai`` + optional ``host.gtf``). Downloaded exactly once;
  shared across all pipeline runs regardless of how many samples or vectors
  are processed. In production (auto mode) this lives at
  ``pipeline_outputs/tla-pip/reference_ncbi/<accession>/``.
* ``<mixed_dir>/<accession>__<plasmid_stem>/`` -- the indexed mixed
  reference (``mixed.fa`` + bwa/faidx index + ``plasmid_contigs.txt``). Built
  once per unique (host, vector) pair; samples sharing a vector reuse it.

When called through ``pipeline.py`` the cache directories are always passed
explicitly and the env-var globals below are not used. The env vars and module-
level defaults are only active when running ``build_reference.py`` standalone.
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import time

from pipeline_io import is_nonempty, get_logger

log = get_logger()

# --- Global cache locations (standalone CLI only; pipeline.py always passes
# explicit host_dir / mixed_dir so these globals are not used in pipeline runs).
HOST_GENOME_DIR = os.environ.get("TLA_HOST_GENOME_DIR", "data/host_genome")
MIXED_REF_DIR = os.environ.get("TLA_MIXED_REF_DIR", "data/mixed_ref")
DEFAULT_HOST_ACCESSION = "GCF_003668045.3"

# bwa index sidecar suffixes that must all be present for a complete index.
_BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")


def _fasta_stem(path):
    """Return the FASTA filename without directory or extension."""
    base = os.path.basename(path)
    for ext in (".fasta", ".fna", ".fa"):
        if base.lower().endswith(ext):
            return base[: -len(ext)]
    return os.path.splitext(base)[0]


def _sanitize(name):
    """Make a string safe to use as a single path component."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")


def fasta_contig_names(fasta_path):
    """Return the list of contig names (first whitespace-delimited token of
    each ``>`` header) in a FASTA file."""
    names = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].split()[0])
    return names


def _index_complete(fasta_path):
    """True if ``fasta_path`` has a complete bwa index + ``.fai``."""
    if not is_nonempty(fasta_path):
        return False
    if not is_nonempty(fasta_path + ".fai"):
        return False
    return all(is_nonempty(fasta_path + suf) for suf in _BWA_INDEX_SUFFIXES)


# NCBI datasets download URL template (informational; logged for provenance).
_NCBI_GENOME_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"


def _download_genome_zip(accession, zip_path, retries=4, backoff=5):
    """Download a genome zip with the NCBI ``datasets`` CLI, retrying on
    transient network failures.

    The ``datasets`` client streams over HTTP/2 and occasionally aborts a
    download mid-transfer (e.g. ``stream error ... INTERNAL_ERROR; received
    from peer``), exiting non-zero even though the accession and arguments are
    valid.  Each attempt writes to a fresh ``zip_path`` (any partial file from
    a failed attempt is removed first) and is retried up to ``retries`` times
    with a linear backoff.

    Raises ``subprocess.CalledProcessError`` if every attempt fails.
    """
    cmd = [
        "datasets", "download", "genome", "accession", accession,
        "--include", "genome,gtf",
        "--filename", zip_path,
    ]
    last_exc = None
    for attempt in range(1, retries + 1):
        if os.path.exists(zip_path):
            os.remove(zip_path)   # discard any partial download before retrying
        try:
            subprocess.check_call(cmd)
            return
        except subprocess.CalledProcessError as exc:
            last_exc = exc
            if attempt < retries:
                wait = backoff * attempt
                log.warning(
                    f"NCBI datasets download failed for {accession} "
                    f"(attempt {attempt}/{retries}, exit {exc.returncode}); "
                    f"retrying in {wait}s. This is usually a transient "
                    f"network/HTTP2 error, not a bad accession."
                )
                time.sleep(wait)
    log.error(
        f"NCBI datasets download failed for {accession} after {retries} "
        f"attempts. Last exit code: {last_exc.returncode}."
    )
    raise last_exc


def ensure_host_genome(accession=DEFAULT_HOST_ACCESSION, host_dir=None, force=False):
    """Download (once) and return the cached host genome FASTA for ``accession``.

    Reuses ``<host_dir>/<accession>/host.fa`` if present (and ``force`` is
    False), otherwise downloads the genome (and GTF, if available) with the
    NCBI ``datasets`` CLI, then runs ``samtools faidx``.

    Parameters
    ----------
    accession : str
        NCBI genome accession (e.g. ``GCF_003668045.3``).
    host_dir : str or None
        Parent directory for the per-accession cache.  ``None`` falls back to
        the module-level ``HOST_GENOME_DIR`` (standalone CLI only).
    force : bool
        When True, delete any existing cache and re-download.

    Returns the path to ``host.fa``.
    """
    host_dir = host_dir if host_dir is not None else HOST_GENOME_DIR
    acc_dir  = os.path.join(host_dir, _sanitize(accession))
    host_fa  = os.path.join(acc_dir, "host.fa")

    if force and os.path.exists(acc_dir):
        log.info(f"Force-download requested; removing existing cache: {acc_dir}")
        shutil.rmtree(acc_dir)

    if is_nonempty(host_fa) and is_nonempty(host_fa + ".fai"):
        log.info(f"Host genome reused from cache: {host_fa}  "
                 f"(accession: {accession})")
        return host_fa

    if is_nonempty(host_fa) and not is_nonempty(host_fa + ".fai"):
        log.warning(f"Host FASTA present but .fai missing; re-indexing: {host_fa}")
        subprocess.check_call(["samtools", "faidx", host_fa])
        return host_fa

    # Detect a stale acc_dir left by a previous interrupted atomic rename
    # (exists but host.fa absent or empty).  Clean it so the download retries.
    if os.path.exists(acc_dir) and not is_nonempty(host_fa):
        log.warning(f"Stale cache directory found (no host.fa); "
                    f"removing and re-downloading: {acc_dir}")
        shutil.rmtree(acc_dir)

    log.info(f"Downloading host genome from NCBI: accession={accession}")
    log.info(f"  source URL (informational): "
             f"{_NCBI_GENOME_URL.format(accession=accession)}")

    # Ensure host_dir (the parent of acc_dir) exists, but do NOT pre-create
    # acc_dir itself -- the atomic rename at the end is the only thing that
    # publishes it, guaranteeing a partial download is never treated as a
    # cache hit.
    os.makedirs(host_dir, exist_ok=True)
    build_dir = acc_dir + ".building"
    if os.path.exists(build_dir):
        log.warning(f"Removing leftover build directory: {build_dir}")
        shutil.rmtree(build_dir)
    os.makedirs(build_dir)

    try:
        zip_path = os.path.join(build_dir, "ncbi_dataset.zip")
        _download_genome_zip(accession, zip_path)
        subprocess.check_call(["unzip", "-o", "-q", zip_path, "-d", build_dir])

        data_root = os.path.join(build_dir, "ncbi_dataset", "data")
        fna_hits = glob.glob(os.path.join(data_root, "**", "*_genomic.fna"), recursive=True)
        if not fna_hits:
            fna_hits = glob.glob(os.path.join(data_root, "**", "*.fna"), recursive=True)
        if not fna_hits:
            raise FileNotFoundError(
                f"No genome FASTA (*.fna) found in NCBI download for {accession}. "
                f"Download zip: {zip_path}")

        tmp_host_fa = os.path.join(build_dir, "host.fa")
        shutil.copyfile(fna_hits[0], tmp_host_fa)

        gtf_hits = glob.glob(os.path.join(data_root, "**", "*.gtf"), recursive=True)
        if gtf_hits:
            shutil.copyfile(gtf_hits[0], os.path.join(build_dir, "host.gtf"))
        else:
            log.warning(f"No GTF in NCBI download for {accession} "
                        f"(continuing -- annotation is optional in V1).")

        subprocess.check_call(["samtools", "faidx", tmp_host_fa])

        # Atomically publish the completed download into the canonical acc_dir.
        if os.path.exists(acc_dir):
            shutil.rmtree(acc_dir)
        os.rename(build_dir, acc_dir)
    except BaseException:
        shutil.rmtree(build_dir, ignore_errors=True)
        raise

    log.info(f"Host genome ready: {host_fa}  (accession: {accession})")
    return host_fa


def build_mixed_reference(goi_fasta, accession=DEFAULT_HOST_ACCESSION,
                          host_dir=None, mixed_dir=None, force=False):
    """Return an indexed mixed (host + GOI) reference, building it if needed.

    The result is cached under ``<mixed_dir>/<accession>__<plasmid_stem>/`` and
    reused on any later call with the same host+GOI pair (so samples sharing a
    vector never re-index).  Pass ``force=True`` to force re-download of the
    host genome before building.

    Returns
    -------
    (mixed_fa_path, plasmid_contigs_path)
    """
    mixed_dir = mixed_dir if mixed_dir is not None else MIXED_REF_DIR
    if not is_nonempty(goi_fasta):
        raise FileNotFoundError(f"GOI/vector FASTA not found or empty: {goi_fasta}")

    host_fa = ensure_host_genome(accession=accession, host_dir=host_dir, force=force)

    # Cache key for this mixed reference: ``<host_accession>__<plasmid_stem>``.
    # The host genome is keyed by accession (the on-disk FASTA is always named
    # ``host.fa``), so the accession -- not the FASTA stem -- identifies the host.
    key = f"{_sanitize(accession)}__{_sanitize(_fasta_stem(goi_fasta))}"
    cache_dir = os.path.join(mixed_dir, key)
    mixed_fa = os.path.join(cache_dir, "mixed.fa")
    plasmid_contigs = os.path.join(cache_dir, "plasmid_contigs.txt")

    # Cache hit: mixed.fa + full index + plasmid_contigs.txt all present.
    if _index_complete(mixed_fa) and is_nonempty(plasmid_contigs):
        log.info(f"Mixed reference already cached: {mixed_fa}")
        return mixed_fa, plasmid_contigs

    # Build the mixed reference: concatenate the host genome and the GOI/vector
    # FASTA, then index. On any failure, wipe the half-built cache dir so the
    # next run rebuilds instead of treating it as a cache hit.
    if os.path.exists(cache_dir):
        shutil.rmtree(cache_dir)
    os.makedirs(cache_dir, exist_ok=True)

    try:
        log.info(f"Building mixed reference [{key}] ...")
        # Concatenate host + GOI; a trailing newline after each source keeps
        # the next header on its own line.
        with open(mixed_fa, "wb") as out:
            for src in (host_fa, goi_fasta):
                with open(src, "rb") as fh:
                    shutil.copyfileobj(fh, out)
                out.write(b"\n")

        log.info("Indexing mixed reference (bwa index) ...")
        subprocess.check_call(["bwa", "index", mixed_fa])
        log.info("Indexing mixed reference (samtools faidx) ...")
        subprocess.check_call(["samtools", "faidx", mixed_fa])

        # Record the GOI/vector contig names so the analysis stage can tag
        # reads as host vs vector.
        with open(plasmid_contigs, "w") as fh:
            for n in fasta_contig_names(goi_fasta):
                fh.write(n + "\n")
    except BaseException:
        shutil.rmtree(cache_dir, ignore_errors=True)
        raise

    log.info(f"Mixed reference ready: {mixed_fa}")
    return mixed_fa, plasmid_contigs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download/reuse the host genome and build a cached mixed "
                    "(host + GOI) reference for the TLA pipeline.")
    parser.add_argument("-g", "--goi_fasta", required=True,
                        help="Vector / GOI reference FASTA to add to the host genome")
    parser.add_argument("-a", "--accession", default=DEFAULT_HOST_ACCESSION,
                        help=f"Host genome NCBI accession (default: {DEFAULT_HOST_ACCESSION})")
    parser.add_argument("--host_dir", default=None,
                        help=f"Host genome cache dir (default: {HOST_GENOME_DIR})")
    parser.add_argument("--mixed_dir", default=None,
                        help=f"Mixed reference cache dir (default: {MIXED_REF_DIR})")
    parser.add_argument("--force", action="store_true",
                        help="Re-download the host genome even if a cache exists")
    args = parser.parse_args()
    mixed_fa, plasmid_contigs = build_mixed_reference(
        args.goi_fasta, accession=args.accession,
        host_dir=args.host_dir, mixed_dir=args.mixed_dir, force=args.force)
    log.info(f"mixed reference: {mixed_fa}")
    log.info(f"plasmid contigs: {plasmid_contigs}")
