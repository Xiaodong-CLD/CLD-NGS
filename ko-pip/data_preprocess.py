'''
This script is a modification .....
'''
import argparse
import subprocess
import os


def fastp_runner(input_fastq_file_1_path, input_fastq_file_2_path):
    trimmed_fastq_file_1_path = input_fastq_file_1_path.replace(".fastq", ".trimmed.fastq")
    trimmed_fastq_file_2_path = input_fastq_file_2_path.replace(".fastq", ".trimmed.fastq")
    fastp_html_report_path = input_fastq_file_1_path.replace(".R1.fastq.gz", ".fastp.html")
    fastp_json_report_path = input_fastq_file_1_path.replace(".R1.fastq.gz", ".fastp.json")


    #Check if fastp has already been executed.
    if os.path.exists(trimmed_fastq_file_1_path) and os.path.exists(trimmed_fastq_file_2_path):
        print("[Takeda Pipeline Message] Read Trimmer: fastp, has already been executed on files {} and {}".format(input_fastq_file_1_path.split('/')[-1], input_fastq_file_2_path.split('/')[-1]))
        return trimmed_fastq_file_1_path, trimmed_fastq_file_2_path

    try:
        print("[Takeda Pipeline Message] Trimming and Quality checking fastq files: {} and {}".format(input_fastq_file_1_path.split('/')[-1], input_fastq_file_2_path.split('/')[-1]))
        bash_output = subprocess.check_output("fastp --in1 " + input_fastq_file_1_path + " --in2 " + input_fastq_file_2_path + " --out1 " + trimmed_fastq_file_1_path + " --out2 " + trimmed_fastq_file_2_path + " -w 16 -j " + fastp_json_report_path + " -h " + fastp_html_report_path, shell=True)

        return trimmed_fastq_file_1_path, trimmed_fastq_file_2_path

    except subprocess.CalledProcessError:
        print("[Takeda Pipeline Message] fastp processing failed.")



def error_correction(trimmed_fastq_file_1_path, trimmed_fastq_file_2_path):
    # Output file names
    corrected_fastq_file_1_path = os.path.splitext(trimmed_fastq_file_1_path)[0] + "_corrected.fastq"
    corrected_fastq_file_2_path = os.path.splitext(trimmed_fastq_file_2_path)[0] + "_corrected.fastq"

    #Check if fastp has already been executed.
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

def data_preprocessing(input_fastq_file_1_path, input_fastq_file_2_path):

    trimmed_fastq_file_1_path, trimmed_fastq_file_2_path = fastp_runner(input_fastq_file_1_path, input_fastq_file_2_path)
    corrected_fastq_file_1_path, corrected_fastq_file_2_path = error_correction(trimmed_fastq_file_1_path, trimmed_fastq_file_2_path)
    
    return corrected_fastq_file_1_path, corrected_fastq_file_2_path

if __name__ == "__main__":
    pass