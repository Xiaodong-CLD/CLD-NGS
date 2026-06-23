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


def fastp_runner(input_fastq_file_1_path, input_fastq_file_2_path,
                 output_dir, naming_prefix, ext=".fastq.gz", report_dir=None):
    os.makedirs(output_dir, exist_ok=True)
    trimmed_fastq_file_1_path = os.path.join(output_dir, f"{naming_prefix}.R1.trimmed{ext}")
    trimmed_fastq_file_2_path = os.path.join(output_dir, f"{naming_prefix}.R2.trimmed{ext}")

    report_dir = report_dir if report_dir is not None else output_dir
    os.makedirs(report_dir, exist_ok=True)
    fastp_html_report_path = os.path.join(report_dir, f"{naming_prefix}.fastp.html")
    fastp_json_report_path = os.path.join(report_dir, f"{naming_prefix}.fastp.json")

    # Check if fastp has already been executed.
    if os.path.exists(trimmed_fastq_file_1_path) and os.path.exists(trimmed_fastq_file_2_path):
        print("[Takeda Pipeline Message] Read Trimmer: fastp, has already been executed on files {} and {}".format(
            os.path.basename(input_fastq_file_1_path), os.path.basename(input_fastq_file_2_path)))
        return trimmed_fastq_file_1_path, trimmed_fastq_file_2_path

    try:
        print("[Takeda Pipeline Message] Trimming and Quality checking fastq files: {} and {}".format(
            os.path.basename(input_fastq_file_1_path), os.path.basename(input_fastq_file_2_path)))
        subprocess.check_call([
            "fastp",
            "--in1", input_fastq_file_1_path,
            "--in2", input_fastq_file_2_path,
            "--out1", trimmed_fastq_file_1_path,
            "--out2", trimmed_fastq_file_2_path,
            "-w", "16",
            "-j", fastp_json_report_path,
            "-h", fastp_html_report_path,
        ])
        return trimmed_fastq_file_1_path, trimmed_fastq_file_2_path

    except subprocess.CalledProcessError:
        print("[Takeda Pipeline Message] fastp processing failed.")
        return False, False



def error_correction(trimmed_fastq_file_1_path, trimmed_fastq_file_2_path):
    # Output file names
    corrected_fastq_file_1_path = os.path.splitext(trimmed_fastq_file_1_path)[0] + "_corrected.fastq"
    corrected_fastq_file_2_path = os.path.splitext(trimmed_fastq_file_2_path)[0] + "_corrected.fastq"

    #Check if bfc has already been executed.
    if os.path.exists(corrected_fastq_file_1_path) and os.path.exists(corrected_fastq_file_2_path):
        print("[Takeda Pipeline Message] Read Correction: bfc, has already been executed on files {} and {}".format(trimmed_fastq_file_1_path.split('/')[-1], trimmed_fastq_file_2_path.split('/')[-1]))
        return corrected_fastq_file_1_path, corrected_fastq_file_2_path
    try:
        # BFC commands for forward and reverse reads
        bfc_cmd_forward = f"bfc {trimmed_fastq_file_1_path} > {corrected_fastq_file_1_path}"
        bfc_cmd_reverse = f"bfc {trimmed_fastq_file_2_path} > {corrected_fastq_file_2_path}"

        # Run BFC on forward reads
        subprocess.run(bfc_cmd_forward, shell=True, check=True)

        # Run BFC on reverse reads
        subprocess.run(bfc_cmd_reverse, shell=True, check=True)

        return corrected_fastq_file_1_path, corrected_fastq_file_2_path
    except subprocess.CalledProcessError:
        print("[Takeda Pipeline Message] bfc processing failed.")


def data_preprocessing(input_fastq_file_1_path, input_fastq_file_2_path,
                       output_dir, naming_prefix, ext=".fastq.gz", report_dir=None):
    trimmed_fastq_file_1_path, trimmed_fastq_file_2_path = fastp_runner(
        input_fastq_file_1_path, input_fastq_file_2_path,
        output_dir, naming_prefix, ext=ext, report_dir=report_dir)
    corrected_fastq_file_1_path, corrected_fastq_file_2_path = error_correction(
        trimmed_fastq_file_1_path, trimmed_fastq_file_2_path)

    return corrected_fastq_file_1_path, corrected_fastq_file_2_path


if __name__ == "__main__":
    pass