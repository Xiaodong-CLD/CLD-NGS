#!/usr/bin/env python

'''
This script is a modification of the script located at:
	https://github.com/ShivamSharma13/genome-assembly-pipeline/blob/master/spades_wrapper.py

The above mentioned script on github is written by me (Shivam Sharma).
'''


import subprocess
import os
import random

def varscan_runner(mpileup_file_path, output_dir=None):
    """Call SNP, micro-indel, and consensus variants with the VarScan mpileup2* family.

    SNP + consensus VCFs land under ``output_dir/snv/`` and the indel VCF under
    ``output_dir/indel/`` (falling back next to the mpileup file when
    ``output_dir`` is None). SNP and indel pass ``--min-var-freq 0.01`` because
    mpileup2* otherwise defaults to 0.20, which drops the low-frequency variants
    the downstream filters still rely on. ``mpileup2indel --output-vcf 1`` emits
    a real, bcftools-normalizable VCF. Returns the consensus VCF path, or False
    on failure.
    """
    base_dir = output_dir if output_dir is not None else os.path.dirname(mpileup_file_path)
    snv_out_dir = os.path.join(base_dir, 'snv')
    indel_out_dir = os.path.join(base_dir, 'indel')
    os.makedirs(snv_out_dir, exist_ok=True)
    os.makedirs(indel_out_dir, exist_ok=True)

    mpileup_basename = os.path.basename(mpileup_file_path)
    output_snp_vcf_file_path = os.path.join(snv_out_dir, mpileup_basename.replace(".mpileup", ".snp.vcf"))
    output_indel_vcf_file_path = os.path.join(indel_out_dir, mpileup_basename.replace(".mpileup", ".indel.vcf"))
    output_consensus_vcf_file_path = os.path.join(snv_out_dir, mpileup_basename.replace(".mpileup", ".consensus.vcf"))

    #Check if varscan has already been executed.
    if os.path.exists(output_consensus_vcf_file_path):
        print("[Takeda Pipeline Message] Varscan has already been executed on sample file {}".format(mpileup_file_path.split('/')[-1]))
        return output_consensus_vcf_file_path

    try:
        # Call SNPs.
        print("[Takeda Pipeline Message] Now running VarScan (SNP) on the mpileup file produced...")
        snp_cmd = "varscan mpileup2snp " + mpileup_file_path + " --min-coverage 1500 --min-var-freq 0.01 > " + output_snp_vcf_file_path
        bash_output = subprocess.check_output(snp_cmd, shell=True)

        # Call micro-indels.
        print("[Takeda Pipeline Message] Now running VarScan (InDel) on the mpileup file produced...")
        print("[Takeda Pipeline Message] VarScan is going to be used for calling micro-InDels (< 50 bases).")
        indel_cmd = "varscan mpileup2indel " + mpileup_file_path + " --min-var-freq 0.01 --output-vcf 1 > " + output_indel_vcf_file_path
        bash_output = subprocess.check_output(indel_cmd, shell=True)

        # Call per-position consensus.
        print("[Takeda Pipeline Message] Now running VarScan for consensus call")
        cns_cmd = "varscan mpileup2cns " + mpileup_file_path + " --min-coverage 1500 > " + output_consensus_vcf_file_path
        bash_output = subprocess.check_output(cns_cmd, shell=True)

        return output_consensus_vcf_file_path

    except subprocess.CalledProcessError:
        print("[Takeda Pipeline Message] VarScan processing failed.")
        return False

if __name__ == "__main__":
	#This script will be called by the main script.
	pass
