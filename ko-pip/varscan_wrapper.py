#!/usr/bin/env python

'''
This script is a modification of the script located at:
	https://github.com/ShivamSharma13/genome-assembly-pipeline/blob/master/spades_wrapper.py

The above mentioned script on github is written by me (Shivam Sharma).
'''


import subprocess
import os
import random

from pipeline_io import is_nonempty

MIN_VAR_FREQ = 0.01  # minimum variant frequency to report in the final VCF

def varscan_runner(mpileup_file_path, snv_dir, indel_dir):
    """Call SNP, micro-indel, and consensus variants with the VarScan mpileup2* family.

    SNP + consensus VCFs land under ``snv_dir`` and the indel VCF under
    ``indel_dir`` (both decided by the orchestrator; in KO-pip these are
    diagnostic intermediates under ``tmp/<sample>/``). All three callers use
    VarScan default for coverage and 0.01 for frequency.
    ``mpileup2indel --output-vcf 1`` emits a real, bcftools-normalizable VCF.
    Returns the consensus VCF path, or False on failure.
    """
    snv_out_dir = snv_dir
    indel_out_dir = indel_dir
    os.makedirs(snv_out_dir, exist_ok=True)
    os.makedirs(indel_out_dir, exist_ok=True)

    mpileup_basename = os.path.basename(mpileup_file_path)
    output_snp_vcf_file_path = os.path.join(snv_out_dir, mpileup_basename.replace(".mpileup", ".snp.vcf"))
    output_indel_vcf_file_path = os.path.join(indel_out_dir, mpileup_basename.replace(".mpileup", ".indel.vcf"))
    output_consensus_vcf_file_path = os.path.join(snv_out_dir, mpileup_basename.replace(".mpileup", ".consensus.vcf"))

    #Check if varscan has already been executed. Require all three VCFs to be
    # non-empty: each is written via a `>` redirect, so an interrupted run can
    # leave a 0-byte file that an existence-only check would wrongly reuse.
    if is_nonempty(output_snp_vcf_file_path) \
            and is_nonempty(output_indel_vcf_file_path) \
            and is_nonempty(output_consensus_vcf_file_path):
        print("[Takeda Pipeline Message] Varscan has already been executed on sample file {}".format(mpileup_file_path.split('/')[-1]))
        return output_consensus_vcf_file_path

    # Temp paths for atomic writes: each VarScan call writes to a `.tmp` file
    # which is renamed into place only after it exits 0, so a killed VarScan
    # never leaves a partial VCF at the canonical path.
    snp_tmp = output_snp_vcf_file_path + ".tmp"
    indel_tmp = output_indel_vcf_file_path + ".tmp"
    consensus_tmp = output_consensus_vcf_file_path + ".tmp"

    # Remove any partial leftovers from a previous failed run so we start clean.
    for p in (output_snp_vcf_file_path, output_indel_vcf_file_path, output_consensus_vcf_file_path,
              snp_tmp, indel_tmp, consensus_tmp):
        if os.path.exists(p):
            os.remove(p)

    try:
        # Call SNPs.
        print("[Takeda Pipeline Message] Now running VarScan (SNP) on the mpileup file produced...")
        snp_cmd = "varscan mpileup2snp " + mpileup_file_path + " --strand-filter 0" + " --min-var-freq " + str(MIN_VAR_FREQ) + " > " + snp_tmp
        bash_output = subprocess.check_output(snp_cmd, shell=True)
        os.replace(snp_tmp, output_snp_vcf_file_path)

        # Call micro-indels.
        print("[Takeda Pipeline Message] Now running VarScan (InDel) on the mpileup file produced...")
        print("[Takeda Pipeline Message] VarScan is going to be used for calling micro-InDels (< 50 bases).")
        indel_cmd = "varscan mpileup2indel " + mpileup_file_path + " --min-var-freq " + str(MIN_VAR_FREQ)    + " --output-vcf 1 > " + indel_tmp
        bash_output = subprocess.check_output(indel_cmd, shell=True)
        os.replace(indel_tmp, output_indel_vcf_file_path)

        # Call per-position consensus.
        print("[Takeda Pipeline Message] Now running VarScan for consensus call")
        cns_cmd = "varscan mpileup2cns " + mpileup_file_path + " > " + consensus_tmp
        bash_output = subprocess.check_output(cns_cmd, shell=True)
        os.replace(consensus_tmp, output_consensus_vcf_file_path)

        return output_consensus_vcf_file_path

    except subprocess.CalledProcessError:
        print("[Takeda Pipeline Message] VarScan processing failed.")
        for p in (output_snp_vcf_file_path, output_indel_vcf_file_path, output_consensus_vcf_file_path,
                  snp_tmp, indel_tmp, consensus_tmp):
            if os.path.exists(p):
                os.remove(p)
        return False

if __name__ == "__main__":
	#This script will be called by the main script.
	pass
