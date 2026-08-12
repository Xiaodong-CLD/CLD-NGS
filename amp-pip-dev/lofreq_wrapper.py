#!/usr/bin/env python

'''
LoFreq indel caller wrapper for the amp-pip-dev pipeline.

Requires a sorted, indexed BAM (produced by samtools_wrapper.py).
Steps:
  1. lofreq indelqual --dindel  — add DINDEL-based indel quality scores to the BAM.
  2. samtools index             — index the indelqual BAM.
  3. lofreq call-parallel       — call indels only.

Output: <sample_dir>/lofreq-output/<sample>.lofreq.indel.vcf
'''

import os
import subprocess


def lofreq_runner(sorted_bam_file_path, reference_file_path):
    sample_dir = "/".join(sorted_bam_file_path.split('/')[:-1])
    lofreq_dir = os.path.join(sample_dir, 'lofreq-output')
    if not os.path.exists(lofreq_dir):
        os.mkdir(lofreq_dir)

    sample_name = sorted_bam_file_path.split('/')[-1].replace(".sorted.bam", "")
    indelqual_bam = os.path.join(lofreq_dir, sample_name + ".indelqual.bam")
    output_vcf = os.path.join(lofreq_dir, sample_name + ".lofreq.indel.vcf")

    if os.path.exists(output_vcf):
        print("[Takeda Pipeline Message] LoFreq output already exists, skipping...")
        return output_vcf

    try:
        print("[Takeda Pipeline Message] Running lofreq indelqual (DINDEL) on {}".format(
            sorted_bam_file_path.split('/')[-1]))
        subprocess.check_output(
            "lofreq indelqual --dindel -f " + reference_file_path
            + " -o " + indelqual_bam + " " + sorted_bam_file_path,
            shell=True,
        )

        subprocess.check_output("samtools index " + indelqual_bam, shell=True)

        print("[Takeda Pipeline Message] Running lofreq call-parallel (indels only) on {}".format(
            indelqual_bam.split('/')[-1]))
        subprocess.check_output(
            "lofreq call-parallel --pp-threads 16 --call-indels --only-indels"
            + " -f " + reference_file_path
            + " -o " + output_vcf
            + " " + indelqual_bam,
            shell=True,
        )
        return output_vcf

    except subprocess.CalledProcessError:
        print("[Takeda Pipeline Message] LoFreq failed on {}.".format(
            sorted_bam_file_path.split('/')[-1]))
        return None


if __name__ == "__main__":
    pass
