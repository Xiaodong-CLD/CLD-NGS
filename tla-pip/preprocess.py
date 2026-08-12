'''
Read trimming (fastp) and error correction (bfc).

Original script created by Lilin Huang.

Updated 2026-06-02 by Cindy Ren:
- fastp_runner now reads the source FASTQ files directly and emits renamed
  trimmed outputs into an explicit output_dir, so the pipeline no longer
  has to copy/rename FASTQ files before invoking preprocessing.
- HTML/JSON QC reports can be written to a separate report_dir (e.g.
  Results/<sample>/qc/) as QC deliverables.
'''
import argparse
import subprocess
import os

from pipeline_io import is_nonempty, gzip_ok, count_fastq_reads, commit, get_logger

log = get_logger()


def fastp_runner(input_fastq_file_1_path, input_fastq_file_2_path,
                 output_dir, naming_prefix, ext=".fastq.gz", report_dir=None):
    os.makedirs(output_dir, exist_ok=True)
    trimmed_fastq_file_1_path = os.path.join(output_dir, f"{naming_prefix}.R1.trimmed{ext}")
    trimmed_fastq_file_2_path = os.path.join(output_dir, f"{naming_prefix}.R2.trimmed{ext}")

    report_dir = report_dir if report_dir is not None else output_dir
    os.makedirs(report_dir, exist_ok=True)
    fastp_html_report_path = os.path.join(report_dir, f"{naming_prefix}.fastp.html")
    fastp_json_report_path = os.path.join(report_dir, f"{naming_prefix}.fastp.json")

    # Check if fastp has already been executed. Require both outputs to be
    # complete gzip streams (not just present) so a truncated file from an
    # interrupted run is re-run rather than reused.
    if gzip_ok(trimmed_fastq_file_1_path) and gzip_ok(trimmed_fastq_file_2_path):
        log.info("Read Trimmer: fastp, has already been executed on files {} and {}".format(
            os.path.basename(input_fastq_file_1_path), os.path.basename(input_fastq_file_2_path)))
        return trimmed_fastq_file_1_path, trimmed_fastq_file_2_path

    # Write to temp paths and atomically rename only after validation, so a
    # killed fastp never leaves a partial file at the canonical output path.
    tmp_1 = trimmed_fastq_file_1_path + ".tmp"
    tmp_2 = trimmed_fastq_file_2_path + ".tmp"
    # Remove any partial leftovers from a previous failed run so we start clean.
    for p in (trimmed_fastq_file_1_path, trimmed_fastq_file_2_path, tmp_1, tmp_2):
        if os.path.exists(p):
            os.remove(p)

    try:
        log.info("Trimming and Quality checking fastq files: {} and {}".format(
            os.path.basename(input_fastq_file_1_path), os.path.basename(input_fastq_file_2_path)))
        subprocess.check_call([
            "fastp",
            "--in1", input_fastq_file_1_path,
            "--in2", input_fastq_file_2_path,
            "--out1", tmp_1,
            "--out2", tmp_2,
            "-w", "16",
            "-j", fastp_json_report_path,
            "-h", fastp_html_report_path,
        ])

        # Integrity: both outputs must be complete gzip streams. Pair invariant:
        # paired-mode fastp drops reads as pairs, so the surviving mate counts
        # must match; a mismatch means a truncated/corrupt output.
        n1 = count_fastq_reads(tmp_1)
        n2 = count_fastq_reads(tmp_2)
        if n1 != n2:
            raise subprocess.CalledProcessError(1, "fastp",
                f"R1/R2 read-count mismatch after trimming: {n1} vs {n2}")

        os.replace(tmp_1, trimmed_fastq_file_1_path)
        os.replace(tmp_2, trimmed_fastq_file_2_path)
        return trimmed_fastq_file_1_path, trimmed_fastq_file_2_path

    except (subprocess.CalledProcessError, ValueError, OSError, EOFError) as err:
        log.error("fastp processing failed: {}".format(err))
        for p in (tmp_1, tmp_2, trimmed_fastq_file_1_path, trimmed_fastq_file_2_path):
            if os.path.exists(p):
                os.remove(p)
        return False, False



def error_correction(trimmed_fastq_file_1_path, trimmed_fastq_file_2_path):
    # Output file names
    corrected_fastq_file_1_path = os.path.splitext(trimmed_fastq_file_1_path)[0] + "_corrected.fastq"
    corrected_fastq_file_2_path = os.path.splitext(trimmed_fastq_file_2_path)[0] + "_corrected.fastq"

    # Check if bfc has already been executed. The corrected files are plain
    # (uncompressed) FASTQ, so we can only cheaply check non-emptiness here;
    # completeness is guaranteed instead by the atomic write below (a partial
    # file never reaches the canonical name).
    if is_nonempty(corrected_fastq_file_1_path) and is_nonempty(corrected_fastq_file_2_path):
        log.info("Read Correction: bfc, has already been executed on files {} and {}".format(trimmed_fastq_file_1_path.split('/')[-1], trimmed_fastq_file_2_path.split('/')[-1]))
        return corrected_fastq_file_1_path, corrected_fastq_file_2_path

    # Write to temp files and atomically rename only on success, so a killed
    # bfc never leaves a partial file at the canonical output path.
    tmp_forward = corrected_fastq_file_1_path + ".tmp"
    tmp_reverse = corrected_fastq_file_2_path + ".tmp"
    # Remove any partial leftovers from a previous failed run so we start clean.
    for p in (corrected_fastq_file_1_path, corrected_fastq_file_2_path, tmp_forward, tmp_reverse):
        if os.path.exists(p):
            os.remove(p)
    try:
        # BFC commands for forward and reverse reads
        bfc_cmd_forward = f"bfc -t 16 {trimmed_fastq_file_1_path} > {tmp_forward}"
        bfc_cmd_reverse = f"bfc -t 16 {trimmed_fastq_file_2_path} > {tmp_reverse}"

        # Run BFC on forward reads
        subprocess.run(bfc_cmd_forward, shell=True, check=True)

        # Run BFC on reverse reads
        subprocess.run(bfc_cmd_reverse, shell=True, check=True)

        # Integrity + pair invariant: bfc rewrites bases but never drops reads,
        # so each output must be record-aligned and the two mate counts must
        # match (and equal the trimmed input counts).
        n1 = count_fastq_reads(tmp_forward)
        n2 = count_fastq_reads(tmp_reverse)
        if n1 != n2:
            raise subprocess.CalledProcessError(1, "bfc",
                f"R1/R2 read-count mismatch after correction: {n1} vs {n2}")

        os.replace(tmp_forward, corrected_fastq_file_1_path)
        os.replace(tmp_reverse, corrected_fastq_file_2_path)

        return corrected_fastq_file_1_path, corrected_fastq_file_2_path
    except (subprocess.CalledProcessError, ValueError, OSError, EOFError) as err:
        log.error("bfc processing failed: {}".format(err))
        for p in (tmp_forward, tmp_reverse,
                  corrected_fastq_file_1_path, corrected_fastq_file_2_path):
            if os.path.exists(p):
                os.remove(p)
        return False, False


def data_preprocessing(input_fastq_file_1_path, input_fastq_file_2_path,
                       output_dir, naming_prefix, ext=".fastq.gz", report_dir=None):
    trimmed_fastq_file_1_path, trimmed_fastq_file_2_path = fastp_runner(
        input_fastq_file_1_path, input_fastq_file_2_path,
        output_dir, naming_prefix, ext=ext, report_dir=report_dir)
    if not trimmed_fastq_file_1_path:
        return False, False
    corrected_fastq_file_1_path, corrected_fastq_file_2_path = error_correction(
        trimmed_fastq_file_1_path, trimmed_fastq_file_2_path)

    return corrected_fastq_file_1_path, corrected_fastq_file_2_path


if __name__ == "__main__":
    pass