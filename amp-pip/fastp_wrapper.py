#!/usr/bin/env python

'''
This script is a modification of the script located at:
	https://github.com/ShivamSharma13/genome-assembly-pipeline/blob/master/spades_wrapper.py

The above mentioned script on github is written by me (Shivam Sharma).
'''



import subprocess
import os
import random

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
        return False, False

if __name__ == "__main__":
	#This script will be called by the main script.
	pass
