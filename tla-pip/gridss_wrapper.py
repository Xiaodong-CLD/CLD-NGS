#!/usr/bin/env python3

"""Run GRIDSS2 on the mixed reference and keep the VCF + assembly BAM.

Runs GRIDSS2 on a sorted BAM aligned to the mixed (host + vector) reference and
keeps both downstream deliverables:

    breakend VCF   (<out_dir>/<sample>_gridss.vcf.gz [+ .tbi])
    assembly BAM   (<out_dir>/<sample>_assembly.bam)

The VCF carries the PASS host<->vector / vector<->vector BND breakends; the
assembly BAM carries the GRIDSS-assembled contigs used to link events and infer
per-site structure. GRIDSS manages its own working directory
(``<out_dir>/*.gridss.working/``) and resumes an interrupted run from it, so this
wrapper skips only when a complete VCF and assembly BAM already exist; otherwise
it lets GRIDSS (re)run, removing only an incomplete final VCF/BAM first (GRIDSS
refuses to overwrite an existing output VCF).

GRIDSS builds its reference cache (``mixed.fa.gridsscache`` / ``.img``) next to
``mixed.fa`` on first use and reuses it across samples. It is heavy
(~15 min/sample with assembly) and scales sub-linearly past ~8 threads. GRIDSS
2.13.2 is tested on Java 8-11; pin ``openjdk=11`` before productionizing.
"""

import argparse
import os
import subprocess

from pipeline_io import bam_ok, gzip_ok, get_logger

log = get_logger()

DEFAULT_THREADS = 8          # GRIDSS scales sub-linearly past ~8 threads
DEFAULT_JVMHEAP = "32g"      # assembly + variant calling heap (per GRIDSS docs)


def _dedup_bam(sorted_bam, out_dir):
    """Write a duplicate-excluded copy of ``sorted_bam`` and return its path.

    GRIDSS extracts breakend evidence from duplicate-flagged reads (it only
    de-duplicates at VF scoring time), so PCR/optical duplicates inflate the
    ``sv.bam`` evidence set and the host<->host background. This drops every
    record with the SAM duplicate flag (0x400) up front via ``samtools view
    -F 0x400`` so GRIDSS never sees them; GRIDSS's own ``ComputeSamTags``
    preprocess step repairs the resulting SA/supplementary inconsistencies.
    Writes ``<out_dir>/<bam_basename>.dedup.bam`` (+ ``.bai``); skips if a valid
    dedup BAM already exists.
    """
    bam_base = os.path.basename(sorted_bam)
    if bam_base.endswith(".bam"):
        bam_base = bam_base[:-len(".bam")]
    dedup_bam = os.path.join(out_dir, f"{bam_base}.dedup.bam")
    if bam_ok(dedup_bam):
        log.info(f"Duplicate-excluded BAM already present, skipping: {dedup_bam}")
        return dedup_bam
    log.info(f"Excluding duplicate-flagged reads (0x400) from {bam_base}.bam ...")
    subprocess.run(
        ["samtools", "view", "-b", "-F", "0x400", "-o", dedup_bam, sorted_bam],
        check=True)
    subprocess.run(["samtools", "index", dedup_bam], check=True)
    if not bam_ok(dedup_bam):
        raise RuntimeError(f"Duplicate exclusion produced a bad BAM: {dedup_bam}")
    log.info(f"Duplicate-excluded BAM -> {dedup_bam}")
    return dedup_bam


def _read_idsv_reads(idsv_metrics_path):
    """Return the ``READS`` count from a GRIDSS ``.idsv_metrics`` file.

    The file is a Picard-style metrics table: a ``## METRICS CLASS`` line, then a
    tab-separated header row, then a single values row. ``READS`` is a named
    column giving the total reads GRIDSS saw. Returns ``None`` if the file is
    missing or the column cannot be located.
    """
    if not os.path.exists(idsv_metrics_path):
        return None
    with open(idsv_metrics_path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    for i, line in enumerate(lines):
        if line.startswith("## METRICS CLASS"):
            header = lines[i + 1].split("\t")
            values = lines[i + 2].split("\t")
            if "READS" in header:
                return int(values[header.index("READS")])
            return None
    return None


def _count_sv_reads(sv_bam_path):
    """Return the read count in the GRIDSS ``.sv.bam`` (extracted SV evidence).

    This BAM holds every read GRIDSS extracted as breakend evidence (split /
    soft-clipped / indel / one-end-anchored / discordant), so its record count is
    the number of reads GRIDSS considered informative. Returns ``None`` if the
    file is missing.
    """
    if not os.path.exists(sv_bam_path):
        return None
    out = subprocess.check_output(["samtools", "view", "-c", sv_bam_path])
    return int(out.strip())


def collect_read_metrics(out_dir, sorted_bam):
    """Return GRIDSS read-accounting metrics for one sample.

    Reads GRIDSS's own working-directory metric files (no re-scan of the full
    input BAM): ``READS`` from ``<bam>.idsv_metrics`` (total reads) and the record
    count of ``<bam>.sv.bam`` (reads GRIDSS extracted as breakend evidence). Both
    live under ``<out_dir>/<bam_basename>.gridss.working/``. Returns
    ``{"total_reads": int|None, "informative_reads": int|None}``.
    """
    bam_base = os.path.basename(sorted_bam)
    working = os.path.join(out_dir, f"{bam_base}.gridss.working")
    return {
        "total_reads": _read_idsv_reads(
            os.path.join(working, f"{bam_base}.idsv_metrics")),
        "informative_reads": _count_sv_reads(
            os.path.join(working, f"{bam_base}.sv.bam")),
    }



def run(sorted_bam, mixed_fa, out_dir, sample_id,
        threads=DEFAULT_THREADS, jvmheap=DEFAULT_JVMHEAP,
        exclude_duplicates=False):
    """Run GRIDSS2 for one sample; return the VCF + assembly BAM paths.

    ``out_dir`` is the GRIDSS working directory (created if needed). Outputs are
    ``<sample_id>_gridss.vcf.gz`` and ``<sample_id>_assembly.bam``. Skips when
    both already validate. When ``exclude_duplicates`` is set, GRIDSS is run on a
    duplicate-excluded copy of ``sorted_bam`` (see ``_dedup_bam``) so PCR/optical
    duplicates do not inflate the SV evidence or host<->host background; read
    metrics then reflect that dedup input. Returns ``{"vcf": <path>, "assembly":
    <path>, "total_reads": int|None, "informative_reads": int|None}`` (the read
    metrics come from GRIDSS's working-dir metric files; see
    ``collect_read_metrics``).
    """
    os.makedirs(out_dir, exist_ok=True)
    vcf = os.path.join(out_dir, f"{sample_id}_gridss.vcf.gz")
    assembly = os.path.join(out_dir, f"{sample_id}_assembly.bam")

    if not bam_ok(sorted_bam):
        raise FileNotFoundError(f"Sorted BAM missing or truncated: {sorted_bam}")
    if not os.path.exists(mixed_fa):
        raise FileNotFoundError(f"Mixed reference FASTA missing: {mixed_fa}")

    gridss_input = _dedup_bam(sorted_bam, out_dir) if exclude_duplicates else sorted_bam

    if gzip_ok(vcf) and bam_ok(assembly):
        log.info(f"GRIDSS outputs already present, skipping: {vcf}")
        return {"vcf": vcf, "assembly": assembly,
                **collect_read_metrics(out_dir, gridss_input)}

    # We only reach here when a complete run is NOT already present. GRIDSS
    # refuses to overwrite an existing output VCF, so remove any final VCF (+
    # index) before (re)running; also drop an incomplete assembly BAM. The
    # *.gridss.working/ subdirs are left in place so GRIDSS resumes.
    for stale in (vcf, vcf + ".tbi"):
        if os.path.exists(stale):
            os.remove(stale)
    if os.path.exists(assembly) and not bam_ok(assembly):
        os.remove(assembly)

    cmd = [
        "gridss",
        "--reference", mixed_fa,
        "--output", vcf,
        "--assembly", assembly,
        "--workingdir", out_dir,
        "--threads", str(threads),
        "--jvmheap", jvmheap,
        gridss_input,
    ]
    log.info(f"Running GRIDSS2 on {os.path.basename(gridss_input)} "
             f"(threads={threads}, jvmheap={jvmheap}) ...")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"GRIDSS run failed for {sample_id}: {exc}") from exc

    if not gzip_ok(vcf):
        raise RuntimeError(f"GRIDSS finished but VCF is missing/truncated: {vcf}")
    if not bam_ok(assembly):
        raise RuntimeError(
            f"GRIDSS finished but assembly BAM is missing/truncated: {assembly}")

    log.info(f"GRIDSS VCF -> {vcf}")
    log.info(f"GRIDSS assembly BAM -> {assembly}")
    metrics = collect_read_metrics(out_dir, gridss_input)
    log.info(f"GRIDSS reads: total={metrics['total_reads']} "
             f"informative={metrics['informative_reads']}")
    return {"vcf": vcf, "assembly": assembly, **metrics}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run GRIDSS2 on a sorted BAM aligned to the mixed "
                    "(host+vector) reference; keeps the breakend VCF + "
                    "assembly BAM.")
    parser.add_argument("-b", "--bam", required=True,
                        help="Sorted BAM aligned to the mixed (host+vector) reference")
    parser.add_argument("-r", "--mixed_fa", required=True,
                        help="Mixed reference FASTA (host + vector)")
    parser.add_argument("-o", "--out_dir", required=True,
                        help="GRIDSS working/output directory (created if needed)")
    parser.add_argument("-n", "--sample_id", required=True,
                        help="Sample name; prefixes the VCF + assembly BAM")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS,
                        help=f"Worker threads (default {DEFAULT_THREADS}; GRIDSS "
                             f"scales sub-linearly past ~8)")
    parser.add_argument("--jvmheap", default=DEFAULT_JVMHEAP,
                        help=f"JVM max heap for assembly + calling (default {DEFAULT_JVMHEAP})")
    parser.add_argument("--exclude-duplicates", action="store_true",
                        help="Drop duplicate-flagged reads (0x400) from the input "
                             "BAM before running GRIDSS, to raise signal-to-noise")
    args = parser.parse_args()
    paths = run(args.bam, args.mixed_fa, args.out_dir, args.sample_id,
                threads=args.threads, jvmheap=args.jvmheap,
                exclude_duplicates=args.exclude_duplicates)
    for name, path in paths.items():
        log.info(f"  {name}: {path}")
