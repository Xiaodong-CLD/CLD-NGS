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

from pipeline_io import is_nonempty, bam_ok, commit

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

    # Skip only if every expected output is present and the BAMs pass a
    # quickcheck (valid header + BGZF EOF block). The non-BAM outputs are
    # written atomically below, so their mere presence implies completeness;
    # mpileup/depth may legitimately be empty (a zero-coverage sample), so we
    # do not require them to be non-empty.
    if bam_ok(sorted_bam_file_path) and os.path.exists(indexed_bam_file_path) \
            and os.path.exists(mpileup_file_path) and os.path.exists(depth_file_path) \
            and bam_ok(collated_bam_file_path):
        print("[Takeda Pipeline Message] Samtools has already been executed on file {}".format(bam_file_path.split('/')[-1]))
        return sorted_bam_file_path, mpileup_file_path, depth_file_path, collated_bam_file_path

    # Temp paths for atomic writes (a killed step never leaves a partial file
    # at the canonical name, so the no-clobber guard above re-runs correctly).
    sorted_tmp = sorted_bam_file_path + ".tmp"
    indexed_tmp = indexed_bam_file_path + ".tmp"
    mpileup_tmp = mpileup_file_path + ".tmp"
    depth_tmp = depth_file_path + ".tmp"
    collated_tmp = collated_bam_file_path + ".tmp"

    # Remove any partial leftovers from a previous failed run so we start clean.
    for p in (sorted_bam_file_path, indexed_bam_file_path, mpileup_file_path,
              depth_file_path, collated_bam_file_path,
              sorted_tmp, indexed_tmp, mpileup_tmp, depth_tmp, collated_tmp):
        if os.path.exists(p):
            os.remove(p)

    try:

        print("[Takeda Pipeline Message] Sorting the .BAM file produced. This sorted file will be used for generating indexed bam (.BAI) file.")
        subprocess.check_call(["samtools", "sort", "-@", "16", "-o", sorted_tmp, bam_file_path])
        if not commit(sorted_tmp, sorted_bam_file_path, validator=bam_ok):
            raise subprocess.CalledProcessError(1, "samtools sort")

        print("[Takeda Pipeline Message] Generating the index for the .BAM file produced. ")
        # Write the index to a temp path, then atomically rename into place.
        subprocess.check_call(["samtools", "index", "-@", "16", sorted_bam_file_path, indexed_tmp])
        if not commit(indexed_tmp, indexed_bam_file_path, validator=is_nonempty):
            raise subprocess.CalledProcessError(1, "samtools index")

        print("[Takeda Pipeline Message] Creating the .mpileup file...")
        # -A  keep anomalous/orphan read pairs (amplicon data).
        # -d 0 disable the per-position depth cap (samtools default is 8000,
        #      which silently truncates this high-depth amplicon to ~1% of
        #      reads and skews every VarScan VarFreq/consensus count).
        # -B  disable BAQ so indel-adjacent bases are not down-weighted, which
        #      preserves indel signal for this KO assay.
        # mpileup may be legitimately empty (zero coverage); a 0 exit means the
        # redirect is complete, so we rename unconditionally on success.
        with open(mpileup_tmp, "w") as fh:
            subprocess.check_call(["samtools", "mpileup", "-f", reference_file_path, "-A", "-d", "0", "-B", sorted_bam_file_path], stdout=fh)
        os.replace(mpileup_tmp, mpileup_file_path)

        print("[Takeda Pipeline Message] Generating coverage files and plots.")
        with open(depth_tmp, "w") as fh:
            subprocess.check_call(["samtools", "depth", sorted_bam_file_path], stdout=fh)
        os.replace(depth_tmp, depth_file_path)

        # Collated (qname-grouped) BAM: mate pairs are adjacent for cheap
        # pair-level iteration. ``collated_tmp + ".collate"`` is samtools'
        # own working-file prefix, distinct from our atomic temp.
        print("[Takeda Pipeline Message] Creating the read-collated .BAM file...")
        subprocess.check_call([
            "samtools", "collate", "-@", "16",
            "-o", collated_tmp,
            sorted_bam_file_path,
            collated_tmp + ".collate",
        ])
        if not commit(collated_tmp, collated_bam_file_path, validator=bam_ok):
            raise subprocess.CalledProcessError(1, "samtools collate")

        return sorted_bam_file_path, mpileup_file_path, depth_file_path, collated_bam_file_path

    except subprocess.CalledProcessError as e:
        print("[Takeda Pipeline Message] SAMtools processing failed: {}".format(e))
        # Clean up partial outputs (final + temp) so the next run doesn't think
        # we succeeded.
        for p in (sorted_bam_file_path, indexed_bam_file_path, mpileup_file_path,
                  depth_file_path, collated_bam_file_path,
                  sorted_tmp, indexed_tmp, mpileup_tmp, depth_tmp, collated_tmp):
            if os.path.exists(p):
                os.remove(p)
        return False, False, False, False

if __name__ == "__main__":
    #This script will be called by the main script.
    pass
