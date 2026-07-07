#!/usr/bin/env python

'''
This script is a modification of the script located at:
    https://github.com/ShivamSharma13/genome-assembly-pipeline/blob/master/spades_wrapper.py

The above mentioned script on github is written by me (Shivam Sharma).
Modified by Tongyu Liu. Remove conversion from sam to bam
'''


import subprocess
import os
import random

def samtools_runner(bam_file_path, reference_file_path):
    sorted_bam_file_path = bam_file_path.replace(".bam", ".sorted.bam")
    indexed_bam_file_path = bam_file_path.replace(".bam", ".sorted.bam.bai")
    mpileup_file_path = bam_file_path.replace(".bam", ".mpileup")
    depth_file_path = bam_file_path.replace(".bam", ".depth.txt")

    #Check if samtools has already been executed.
    if os.path.exists(mpileup_file_path):
        print("[Takeda Pipeline Message] Samtools has already been executed on file {}".format(bam_file_path.split('/')[-1]))
        return bam_file_path, sorted_bam_file_path, mpileup_file_path

    try:

        print("[Takeda Pipeline Message] Sorting the .BAM file produced. This sorted file will be used for generating indexed bam (.BAI) file.")
        bash_output = subprocess.check_output("samtools sort -@ 16 -O bam -o " + sorted_bam_file_path + " " + bam_file_path, shell = True)

        print("[Takeda Pipeline Message] Generating the index for the .BAM file produced. ")
        bash_output = subprocess.check_output("samtools index -@ 16 -b " + sorted_bam_file_path + " > " + indexed_bam_file_path, shell = True)

        print("[Takeda Pipeline Message] Creating the .mpileup file...")
        # -d 0: uncapped max depth for deep amplicon-seq assay.
        # See experiments/2026-06-mpileup-depth-cap/. Do NOT copy to WGS/panels.
        bash_output = subprocess.check_output("samtools mpileup -d 0 -f " + reference_file_path + " " + sorted_bam_file_path + " -A > " + mpileup_file_path, shell = True)

        print("[Takeda Pipeline Message] Generating coverage files and plots.")
        bash_output = subprocess.check_output("samtools depth -a " +  sorted_bam_file_path + " > " + depth_file_path, shell = True)

        return bam_file_path, sorted_bam_file_path, mpileup_file_path

    except subprocess.CalledProcessError:
        print("[Takeda Pipeline Message] SAMtools processing failed.")
        return False, False, False

if __name__ == "__main__":
    #This script will be called by the main script.
    pass
