#!/usr/bin/env python

'''
This script is a modification of the script located at:
	https://github.com/ShivamSharma13/genome-assembly-pipeline/blob/master/spades_wrapper.py

The above mentioned script on github is written by me (Shivam Sharma).
'''


import subprocess
import os
import random

def varscan_runner(mpileup_file_path):
    #Create a directory for varscan results.
    if not os.path.exists(os.path.join("/".join(mpileup_file_path.split('/')[:-1]), 'varscan-output')):
        os.mkdir(os.path.join("/".join(mpileup_file_path.split('/')[:-1]), 'varscan-output'))

    output_snp_vcf_file_path = os.path.join("/".join(mpileup_file_path.split('/')[:-1]), 'varscan-output', mpileup_file_path.split('/')[-1].replace(".mpileup", ".snp.vcf"))
    output_indel_vcf_file_path = os.path.join("/".join(mpileup_file_path.split('/')[:-1]), 'varscan-output', mpileup_file_path.split('/')[-1].replace(".mpileup", ".indel.vcf"))
    output_consensus_vcf_file_path = os.path.join("/".join(mpileup_file_path.split('/')[:-1]), 'varscan-output', mpileup_file_path.split('/')[-1].replace(".mpileup", ".consensus.vcf"))

    #Check if varscan has already been executed.
    if os.path.exists(output_consensus_vcf_file_path):
        print("[Takeda Pipeline Message] Varscan has already been executed on sample file {}".format(mpileup_file_path.split('/')[-1]))
        return output_consensus_vcf_file_path

    try:
        print("[Takeda Pipeline Message] Now running VarScan (SNP) on the PILEUP file produced...")
        bash_output = subprocess.check_output("varscan pileup2snp " + mpileup_file_path + " --min-coverage 1500 > " + output_snp_vcf_file_path, shell = True)

        print("[Takeda Pipeline Message] Now running VarScan (InDel) on the PILEUP file produced...")
        print("[Takeda Pipeline Message] VarScan is going to be used for calling micro-InDels (< 50 bases).")
        bash_output = subprocess.check_output("varscan pileup2indel " + mpileup_file_path + " > " + output_indel_vcf_file_path, shell = True)

        print("[Takeda Pipeline Message] Now running VarScan for consensus call")
        bash_output = subprocess.check_output("varscan pileup2cns " + mpileup_file_path + " > " + output_consensus_vcf_file_path, shell = True)

        return output_consensus_vcf_file_path

    except subprocess.CalledProcessError:
        print("[Takeda Pipeline Message] VarScan processing failed.")
        return False

if __name__ == "__main__":
	#This script will be called by the main script.
	pass
