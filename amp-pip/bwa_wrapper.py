#!/usr/bin/env python

'''
This script is a modification of the script located at:
    https://github.com/ShivamSharma13/genome-assembly-pipeline/blob/master/spades_wrapper.py

The above mentioned script on github is written by me (Shivam Sharma).
Modified by Tongyu Liu, change bwa aln to bwa mem
'''


import subprocess
import os
import random

def bwa_runner(reference_file_path, fastq_file_1_path, fastq_file_2_path, disable_trimming):

    if disable_trimming:
        final_output_file = fastq_file_1_path.replace(".R1.fastq.gz", ".bam")
    else:
        final_output_file = fastq_file_1_path.replace(".R1.trimmed.fastq.gz", ".bam")
    temp_output_file = fastq_file_1_path.replace(".R1.trimmed.fastq.gz", ".trimmed.fq.gz.tmp")

    #Check if BWA has already been executed.
    if os.path.exists(final_output_file):
        print("[Takeda Pipeline Message] BWA has already run on sample {}, skippig bwa alignment process...".format(fastq_file_1_path.split('/')[-1]))
        return final_output_file

    try:
        print("[Takeda Pipeline Message] Running BWA index and alignment...")
        bash_output = subprocess.check_output(["bwa", "index", reference_file_path])

        print("[Takeda Pipeline Message] Running BWA alignment...")
        bash_output = subprocess.check_output("bwa mem -t 6 " + reference_file_path + " " + fastq_file_1_path + " " +  fastq_file_2_path + " | samtools sort -@ 4 -O bam -T " + temp_output_file + " -o " + final_output_file + " -"
, shell=True)

        return final_output_file

    except subprocess.CalledProcessError:
        print("[Takeda Pipeline Message] BWA failed.")
        return False

if __name__ == "__main__":
    #This script will be called by the main script.
    pass
