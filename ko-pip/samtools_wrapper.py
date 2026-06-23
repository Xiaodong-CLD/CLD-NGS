#!/usr/bin/env python

'''
This script is a modification of the script located at:
    https://github.com/ShivamSharma13/genome-assembly-pipeline/blob/master/spades_wrapper.py

The above mentioned script on github is written by me (Shivam Sharma).
Modified by Tongyu Liu. Remove conversion from sam to bam

Modified by Cindy Ren. 
Corrected samtools commands with updated flags and added error handling.
'''


import subprocess
import os
import random

def samtools_runner(bam_file_path, reference_file_path):
    """Produce coord-sorted + indexed BAM, mpileup, depth, and a qname-grouped
    (samtools collate) BAM for pair-level evidence classification.

    Returns
    -------
    (sorted_bam, mpileup, depth, collated_bam)
    """
    sorted_bam_file_path = bam_file_path.replace(".bam", ".sorted.bam")
    indexed_bam_file_path = bam_file_path.replace(".bam", ".sorted.bam.bai")
    mpileup_file_path = bam_file_path.replace(".bam", ".mpileup")
    depth_file_path = bam_file_path.replace(".bam", ".depth.txt")
    collated_bam_file_path = bam_file_path.replace(".bam", ".collated.bam")

    # Skip only if every expected output exists and is non-empty.
    # Previously this checked just the mpileup, which could be left as a 0-byte
    # file by a failed `>` redirect and incorrectly mark the step "done".
    def _ok(p):
        return os.path.exists(p) and os.path.getsize(p) > 0

    if _ok(sorted_bam_file_path) and _ok(indexed_bam_file_path) \
            and _ok(mpileup_file_path) and _ok(depth_file_path) \
            and _ok(collated_bam_file_path):
        print("[Takeda Pipeline Message] Samtools has already been executed on file {}".format(bam_file_path.split('/')[-1]))
        return sorted_bam_file_path, mpileup_file_path, depth_file_path, collated_bam_file_path

    # Remove any partial leftovers from a previous failed run so we start clean.
    for p in (sorted_bam_file_path, indexed_bam_file_path, mpileup_file_path,
              depth_file_path, collated_bam_file_path):
        if os.path.exists(p):
            os.remove(p)

    try:

        print("[Takeda Pipeline Message] Sorting the .BAM file produced. This sorted file will be used for generating indexed bam (.BAI) file.")
        subprocess.check_call(["samtools", "sort", "-@", "16", "-o", sorted_bam_file_path, bam_file_path])

        print("[Takeda Pipeline Message] Generating the index for the .BAM file produced. ")
        subprocess.check_call(["samtools", "index", "-@", "16", sorted_bam_file_path])

        print("[Takeda Pipeline Message] Creating the .mpileup file...")
        # -A  keep anomalous/orphan read pairs (amplicon data).
        # -d 0 disable the per-position depth cap (samtools default is 8000,
        #      which silently truncates this high-depth amplicon to ~1% of
        #      reads and skews every VarScan VarFreq/consensus count).
        # -B  disable BAQ so indel-adjacent bases are not down-weighted, which
        #      preserves indel signal for this KO assay.
        with open(mpileup_file_path, "w") as fh:
            subprocess.check_call(["samtools", "mpileup", "-f", reference_file_path, "-A", "-d", "0", "-B", sorted_bam_file_path], stdout=fh)

        print("[Takeda Pipeline Message] Generating coverage files and plots.")
        with open(depth_file_path, "w") as fh:
            subprocess.check_call(["samtools", "depth", sorted_bam_file_path], stdout=fh)

        # Collated (qname-grouped) BAM: mate pairs are adjacent for cheap
        # pair-level iteration.
        print("[Takeda Pipeline Message] Creating the read-collated .BAM file...")
        subprocess.check_call([
            "samtools", "collate", "-@", "16",
            "-o", collated_bam_file_path,
            sorted_bam_file_path,
            collated_bam_file_path + ".tmp",
        ])

        return sorted_bam_file_path, mpileup_file_path, depth_file_path, collated_bam_file_path

    except subprocess.CalledProcessError as e:
        print("[Takeda Pipeline Message] SAMtools processing failed: {}".format(e))
        # Clean up partial outputs so the next run doesn't think we succeeded.
        for p in (sorted_bam_file_path, indexed_bam_file_path, mpileup_file_path,
                  depth_file_path, collated_bam_file_path):
            if os.path.exists(p):
                os.remove(p)
        return False, False, False, False

if __name__ == "__main__":
    #This script will be called by the main script.
    pass
