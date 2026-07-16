"""
Alignment + assembly helpers.

Created by Lilin Huang.
Updated by Cindy Ren 2026-06-02.

Changes from previous version version:
- `bwa mem` is now invoked with `-R` to embed a per-sample read group (@RG)
  directly into the SAM/BAM output. This removes the need for a downstream
  Picard `AddOrReplaceReadGroups` step (and the extra samtools sort/index
  that followed it) in `gatk_variant_calling.sh`, simplifying the workflow.
"""

import subprocess
import os
from samtools_wrapper import samtools_runner
from pipeline_io import bam_ok, commit

def alignment(processed_fastq_file_1_path, processed_fastq_file_2_path,reference_path,sample_dir):
    sample_name = os.path.basename(processed_fastq_file_1_path).split(".")[0]
    sam_file_path = os.path.join(sample_dir,f'{sample_name}.sam')
    bam_file_path = os.path.join(sample_dir,f'{sample_name}.bam')
    bam_tmp_path = bam_file_path + ".tmp"

    # Check if the output file already exists, if so skip the alignment.
    # Require a complete BAM (valid header + BGZF EOF block): a truncated BAM
    # from an interrupted `samtools view` must be re-run, not reused.
    if bam_ok(bam_file_path):
        print("[Takeda Pipeline Message] Alignment has already been executed. Skip the alignment.")
        return bam_file_path

    # Remove any partial leftovers from a previous failed run so we start clean.
    for p in (sam_file_path, bam_file_path, bam_tmp_path):
        if os.path.exists(p):
            os.remove(p)
    try:
        # Run alignment. Embed a read group (@RG) directly so downstream GATK
        # tools have the required RG tag without a separate Picard step.
        read_group = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tLB:{sample_name}\\tPL:illumina\\tPU:unit1"
        bwa_command = ["bwa", "mem", "-R", read_group, reference_path, processed_fastq_file_1_path, processed_fastq_file_2_path]
        # Open sam output file
        with open(sam_file_path,"w") as output_file:
            # Run BWA
            subprocess.run(bwa_command, stdout=output_file, check = True)

        # Convert sam file to bam file needed for later analysis. Write to a
        # temp path and atomically rename only after a quickcheck passes.
        samtools_command = ["samtools","view","-S","-b",sam_file_path]
        # Open bam output file
        with open(bam_tmp_path, "w") as output_file:
            # Run samtools
            subprocess.run(samtools_command, stdout=output_file, check=True)

        if not commit(bam_tmp_path, bam_file_path, validator=bam_ok):
            raise subprocess.CalledProcessError(1, "samtools view")
        return bam_file_path

    except subprocess.CalledProcessError:
        print("[Takeda Pipeline Message] Alignment failed.")
        for p in (sam_file_path, bam_file_path, bam_tmp_path):
            if os.path.exists(p):
                os.remove(p)
        return False



def assembly(processed_fastq_file_1_path, processed_fastq_file_2_path,reference_path,sample_dir):
    
    bam_file_path = alignment(processed_fastq_file_1_path, processed_fastq_file_2_path,reference_path,sample_dir)
    sorted_bam_file_path, mpileup_file_path, depth_file_path, collated_bam_file_path = samtools_runner(bam_file_path, reference_path)
    
    return sorted_bam_file_path, mpileup_file_path, depth_file_path, collated_bam_file_path
    
if __name__ == "__main__":
    pass