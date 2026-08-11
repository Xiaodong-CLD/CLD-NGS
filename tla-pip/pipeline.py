#!/usr/bin/env python3

"""TLA-seq pipeline (V1) -- pre-analysis orchestration.

For each sample in the xlsx sample sheet this wires together the independent
modules:

    build mixed reference (host download + vector, cached)
        -> preprocess (fastp trimming + BFC error correction)
        -> alignment (bwa mem vs mixed reference, sorted + indexed BAM)
        -> host + vector coverage
        -> GRIDSS2 integration-site calling
        -> vector variant calling
        -> aggregate run report

Input directory layout (mirrors the ko-pip convention):

    input_dir/
      Fastq/                 paired .fastq[.gz] / .fq[.gz]
      sample_xlsx/*.xlsx     sample sheet (columns below)

Sample sheet columns (required):
    Sample ID (unique), fastq (forward), fastq (reverse)
"""

import argparse
import logging
import os
import pathlib
import sys

import pandas as pd

from build_reference import (build_mixed_reference, ensure_host_genome,
                             DEFAULT_HOST_ACCESSION)
from preprocess import data_preprocessing
from align import alignment, samtools_runner
import host_coverage
import gridss_wrapper
import gridss_sites
import gridss_split_report
import discovery_sites
import vector_variants
from aggregate_run_report import aggregate_run_report
from pipeline_io import get_logger, bam_ok, is_nonempty

log = get_logger()


def _add_file_log_handler(log_path):
    """Attach a file handler to the pipeline logger (idempotent)."""
    logger = logging.getLogger("tla_pip")
    for h in logger.handlers:
        if isinstance(h, logging.FileHandler) and h.baseFilename == os.path.abspath(log_path):
            return  # already attached
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    fh = logging.FileHandler(log_path)
    fh.setFormatter(logging.Formatter(
        "[%(asctime)s] [Takeda Pipeline Message] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"))
    logger.addHandler(fh)
    log.info(f"Pipeline log: {log_path}")


def _fastq_ext(fastq_name):
    """Return the source FASTQ extension (.fastq[.gz] / .fq[.gz]) so trimmed
    outputs keep the same compression."""
    if ".fastq" in fastq_name:
        base = ".fastq"
    else:
        base = ".fq"
    suffix = fastq_name.split(base)[1]  # '' or '.gz'
    return f"{base}{suffix}"


def _write_gridss_summary(path, sample_id, metrics, bp_counts, split_counts):
    """Write a one-row TSV summary for one sample.

    Consolidates GRIDSS read metrics, per-category breakpoint counts, NlaIII
    restriction-site statistics, and resolved/unresolved junction counts.
    Written atomically via a .tmp file.
    """
    def _pct(num, denom):
        return round(100.0 * num / denom, 2) if denom else 0.0

    n_hv  = bp_counts.get("integration_sites", 0) or 0
    n_hv_both = bp_counts.get("n_host_vector_both_ends_nlaiii", 0) or 0
    n_hh  = bp_counts.get("host_host_bnd", 0) or 0
    n_hh_both = bp_counts.get("n_host_host_both_ends_nlaiii", 0) or 0
    n_vv  = bp_counts.get("vector_sv", 0) or 0

    columns = [
        "sample_id", "n_total_reads", "n_sa_tagged_reads",
        "n_host_vector_breakpoints", "n_host_vector_both_ends_nlaiii",
        "pct_host_vector_both_ends_nlaiii", "n_host_vector_breakpoints_retained",
        "n_host_host_breakpoints", "n_host_host_both_ends_nlaiii",
        "pct_host_host_both_ends_nlaiii",
        "n_vector_vector_breakpoints", "n_vector_vector_both_ends_nlaiii",
        "pct_vector_vector_both_ends_nlaiii",
        "n_resolved_sites", "n_unresolved_junctions",
    ]
    values = [
        sample_id,
        metrics.get("total_reads", ""),
        metrics.get("informative_reads", ""),
        n_hv,
        n_hv_both,
        _pct(n_hv_both, n_hv),
        n_hv - n_hv_both,
        n_hh,
        n_hh_both,
        _pct(n_hh_both, n_hh),
        n_vv,
        0,   # vector_sv rows are not NlaIII-annotated
        0.0,
        split_counts.get("n_resolved_sites", ""),
        split_counts.get("n_unresolved_junctions", ""),
    ]
    tmp_path = path + ".tmp"
    with open(tmp_path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        fh.write("\t".join(str(v) for v in values) + "\n")
    os.replace(tmp_path, path)
    return path


def _resolve_reference_dirs(references_dir, output_dir, accession):
    """Return (host_dir, mixed_dir), creating them as needed.

    Local mode  (``references_dir`` supplied): the path must already exist and
    is used as the shared cache root (development / reproducible testing).
    Host genome is expected under ``host_genome/``; mixed refs under
    ``mixed_ref/``.

    Auto mode  (``references_dir`` is None): the host genome is cached in a
    shared directory next to the run output directory so it is downloaded only
    once per accession regardless of how many runs are executed.
    Cache path: ``<parent-of-output_dir>/reference_ncbi/<accession>/``.
    The mixed reference (vector-specific) remains per-run under
    ``<output_dir>/reference/mixed_reference/``.
    """
    if references_dir is not None:
        if not os.path.isdir(references_dir):
            raise ValueError(
                f"--references-dir does not exist or is not a directory: "
                f"{references_dir}")
        host_dir  = os.path.join(references_dir, "host_genome")
        mixed_dir = os.path.join(references_dir, "mixed_ref")
        log.info(f"Reference mode: local (development/testing override)")
        log.info(f"  host_dir : {host_dir}")
        log.info(f"  mixed_dir: {mixed_dir}")
    else:
        run_parent = os.path.dirname(os.path.abspath(output_dir))
        host_dir  = os.path.join(run_parent, "reference_ncbi")
        mixed_dir = os.path.join(output_dir, "reference", "mixed_reference")
        log.info(f"Reference mode: shared NCBI cache + per-run mixed reference")
        log.info(f"  host_dir : {host_dir}/{accession}/  (accession: {accession})")
        log.info(f"  mixed_dir: {mixed_dir}")
    # host_dir is created by ensure_host_genome; only mixed_dir needs mkdir here.
    os.makedirs(mixed_dir, exist_ok=True)
    return host_dir, mixed_dir


def run(sample_sheet, fastq_dir, vector_fasta, output_dir,
        accession=DEFAULT_HOST_ACCESSION, references_dir=None,
        force_download=False):
    # Read sample sheet (.tsv/.txt or .xlsx).
    if sample_sheet.lower().endswith((".tsv", ".txt")):
        sheet = pd.read_csv(sample_sheet, sep="\t").fillna("")
    else:
        sheet = pd.read_excel(sample_sheet).fillna("")

    required_cols = ["Sample ID (unique)", "fastq (forward)", "fastq (reverse)"]
    missing = [c for c in required_cols if c not in sheet.columns]
    if missing:
        raise ValueError(f"Sample sheet is missing required column(s): {missing}")

    host_dir, mixed_dir = _resolve_reference_dirs(references_dir, output_dir, accession)

    # Route pipeline log into the output directory.
    _add_file_log_handler(os.path.join(output_dir, "logs", "pipeline.log"))

    tmp_root = os.path.join(output_dir, "tmp")
    # Samples written directly under output_dir/<sample_id>/.
    results_root = output_dir
    pathlib.Path(tmp_root).mkdir(parents=True, exist_ok=True)
    pathlib.Path(results_root).mkdir(parents=True, exist_ok=True)

    sample_ids = []
    failed_ids: set[str] = set()
    for _, row in sheet.iterrows():
        sample_id = str(row["Sample ID (unique)"]).strip()
        fq1_name = str(row["fastq (forward)"]).strip()
        fq2_name = str(row["fastq (reverse)"]).strip()
        if not sample_id:
            log.warning("Skipping a row with no Sample ID.")
            continue

        log.info(f"==== Sample: {sample_id} ====")

        fq1_path = os.path.join(fastq_dir, fq1_name)
        fq2_path = os.path.join(fastq_dir, fq2_name)

        # 1. Mixed reference (host download + vector), cached under references_dir.
        mixed_fa, plasmid_contigs = build_mixed_reference(
            vector_fasta, accession=accession,
            host_dir=host_dir, mixed_dir=mixed_dir, force=force_download)

        # Per-sample result + tmp dirs.
        sample_results = os.path.join(results_root, sample_id)
        qc_dir = os.path.join(sample_results, "qc")
        align_dir = os.path.join(sample_results, "alignment")
        tmp_sample = os.path.join(tmp_root, sample_id)
        for d in (qc_dir, align_dir, tmp_sample):
            pathlib.Path(d).mkdir(parents=True, exist_ok=True)

        # 2. Preprocess: fastp trimming + BFC error correction. The
        #    BFC-corrected reads are what feed alignment.
        #    Skip guard: if the sorted BAM and fastp HTML are already
        #    present and complete, preprocessing + alignment are
        #    unnecessary even if the tmp FASTQs were cleaned up.
        ext = _fastq_ext(fq1_name)
        flagstat_out = os.path.join(align_dir, f"{sample_id}.flagstats.tsv")
        _expected_sorted_bam = os.path.join(tmp_sample, f"{sample_id}.sorted.bam")
        _fastp_html           = os.path.join(qc_dir, f"{sample_id}.fastp.html")
        if bam_ok(_expected_sorted_bam) and is_nonempty(_fastp_html):
            log.info(f"{sample_id}: sorted BAM and fastp HTML already present; "
                     "skipping preprocessing and alignment.")
            sorted_bam = _expected_sorted_bam
            flagstat   = flagstat_out
        else:
            corrected_1, corrected_2 = data_preprocessing(
                fq1_path, fq2_path,
                output_dir=tmp_sample, naming_prefix=sample_id,
                ext=ext, report_dir=qc_dir)
            if not corrected_1:
                log.error(f"Preprocessing failed for {sample_id}; skipping sample.")
                failed_ids.add(sample_id)
                continue

            # 3. Alignment vs the mixed reference, then sort/index + flagstat.
            bam = alignment(corrected_1, corrected_2, mixed_fa, tmp_sample)
            if not bam:
                log.error(f"Alignment failed for {sample_id}; skipping sample.")
                failed_ids.add(sample_id)
                continue
            # Alignment stats (TSV) are persistent deliverables
            # under Results/<sample>/alignment/, prefixed with the sample name.
            sorted_bam, flagstat = samtools_runner(
                bam, flagstat_path=flagstat_out)
            if not sorted_bam:
                log.error(f"BAM sort/index failed for {sample_id}; skipping sample.")
                failed_ids.add(sample_id)
                continue

        log.info(f"{sample_id}: sorted BAM -> {sorted_bam}")
        log.info(f"{sample_id}: plasmid contigs -> {plasmid_contigs}")
        log.info(f"{sample_id}: alignment stats -> {flagstat}")

        # 4. Host coverage (bedGraph + lane plot + unplaced-scaffold diagnostic).
        #    Vector characterization (VarScan + GRIDSS overlay plot) runs after
        #    GRIDSS so the combined breakpoints file is available.
        host_fai = ensure_host_genome(accession=accession, host_dir=host_dir, force=False) + ".fai"
        coverage_dir = os.path.join(sample_results, "coverage")
        pathlib.Path(coverage_dir).mkdir(parents=True, exist_ok=True)
        cov_prefix = os.path.join(coverage_dir, sample_id)
        try:
            cov_res = host_coverage.run(
                sorted_bam, host_fai, plasmid_contigs, cov_prefix)
        except Exception as exc:
            log.error(f"Coverage plot generation failed for {sample_id}: {exc}")
            failed_ids.add(sample_id)
            continue

        # 5. Integration-site identification (GRIDSS2):
        #      - 4a gridss_wrapper: GRIDSS run on the mixed ref (VCF + assembly
        #        BAM) into the sample tmp dir; also collects read metrics.
        #      - 4b gridss_sites: classify PASS breakends into the four tables
        #        under Results/<sample>/integration/.
        #    A combined run summary (4a read metrics + 4b junction counts) is
        #    written to Results/<sample>/integration/<sample>.gridss_summary.tsv.
        integration_dir = os.path.join(sample_results, "integration")
        gridss_out = os.path.join(tmp_sample, "gridss")
        pathlib.Path(integration_dir).mkdir(parents=True, exist_ok=True)
        breakpoints_path = None
        try:
            gridss_res = gridss_wrapper.run(
                sorted_bam, mixed_fa, gridss_out, sample_id)
            sites_res = gridss_sites.run(
                gridss_res["vcf"], plasmid_contigs,
                os.path.join(integration_dir, sample_id),
                mixed_fa=mixed_fa)
            split_res = gridss_split_report.run(
                sites_res["outputs"]["breakpoints"],
                os.path.join(integration_dir, sample_id),
            )
            disc_counts = discovery_sites.run(
                cov_res["bedgraph"],
                split_res["resolved"],
                split_res["unresolved"],
                os.path.join(integration_dir, sample_id),
            )
            _write_gridss_summary(
                os.path.join(integration_dir, f"{sample_id}.summary.tsv"),
                sample_id, gridss_res, sites_res["counts"], split_res["counts"])
            breakpoints_path = sites_res["outputs"]["breakpoints"]
        except Exception as exc:
            log.error(f"Integration-site identification failed for {sample_id}: {exc}")
            failed_ids.add(sample_id)
            continue

        # 5b. Within-vector characterization (VarScan + GRIDSS vector-vector
        #     breakpoints + per-base coverage plot with overlay).
        try:
            vec_char_dir = os.path.join(sample_results, "vector_characterization")
            pathlib.Path(vec_char_dir).mkdir(parents=True, exist_ok=True)
            vector_variants.run(
                sorted_bam, mixed_fa, plasmid_contigs,
                os.path.join(vec_char_dir, sample_id),
                gridss_breakpoints=breakpoints_path,
                tmp_prefix=os.path.join(tmp_sample, "vector_characterization", sample_id),
                coverage_dir=coverage_dir,
            )
        except Exception as exc:
            log.error(f"Vector characterization failed for {sample_id}: {exc}")
            failed_ids.add(sample_id)

        sample_ids.append(sample_id)

    if failed_ids:
        log.warning(f"TLA pipeline finished with {len(failed_ids)} failed sample(s): {sorted(failed_ids)}")
    log.info("TLA pre-analysis pipeline complete.")

    # Aggregate run-level report across all processed samples.
    final_report_dir = os.path.join(output_dir, "final_report")
    aggregate_run_report(results_root, final_report_dir, sorted(sample_ids))
    return len(failed_ids)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TLA-seq pipeline (V1).")
    parser.add_argument("-s", "--sample-sheet", required=True,
                        help="Sample sheet (.tsv or .xlsx)")
    parser.add_argument("-f", "--fastq-dir", required=True,
                    help="Directory containing the FASTQ files listed in the sample sheet")
    parser.add_argument("-v", "--vector-fasta", required=True,
                        help="Vector/GOI FASTA (shared across all samples in this run)")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="Output directory; one subfolder per sample is created here")
    parser.add_argument("-r", "--references-dir", default=None,
                        help="Local host-genome and mixed-ref cache directory "
                             "(development/testing override). Must already exist. "
                             "When omitted, the host genome is cached under "
                             "<parent-of-output-dir>/reference_ncbi/<accession>/ "
                             "(shared across runs) and the mixed reference is "
                             "placed per-run under <output-dir>/reference/mixed_reference/.")
    parser.add_argument("-a", "--accession", default=DEFAULT_HOST_ACCESSION,
                        help=f"Host genome NCBI accession (default: {DEFAULT_HOST_ACCESSION})")
    parser.add_argument("--force-download", action="store_true",
                        help="Re-download the host genome even if a cached copy exists")
    args = parser.parse_args()

    n_failed = run(
        sample_sheet=args.sample_sheet,
        fastq_dir=args.fastq_dir,
        vector_fasta=args.vector_fasta,
        output_dir=args.output_dir,
        accession=args.accession,
        references_dir=args.references_dir,
        force_download=args.force_download,
    )
    sys.exit(1 if n_failed else 0)
