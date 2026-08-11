"""
Alignment helpers for the TLA pipeline.

Created by Lilin Huang.
Adapted for TLA by Cindy Ren.

- `bwa mem` aligns reads against the mixed (host + vector) reference. Default
  bwa-mem behaviour is preserved so supplementary alignments, soft clips,
  MAPQ, CIGAR and SA tags survive into the BAM -- these are the evidence the
  TLA analysis stage relies on.
- `samtools_runner` (previously a separate `samtools_wrapper` module) sorts +
  marks duplicates + indexes the BAM and writes a flagstat summary. PCR/optical
  duplicates are flagged (0x400) via `fixmate -m` + `markdup` so downstream
  stages can skip them (TLA libraries are PCR-heavy and duplicates pile up in
  background coverage, which would otherwise inflate junction-fraction
  denominators). A genome-wide mpileup would be huge against a whole host genome
  and is not produced; per-region depth/coverage is computed later in the
  analysis stage.
"""

import argparse
import subprocess
import os
from pipeline_io import is_nonempty, bam_ok, commit, get_logger

log = get_logger()

def alignment(processed_fastq_file_1_path, processed_fastq_file_2_path,reference_path,sample_dir):
    sample_name = os.path.basename(processed_fastq_file_1_path).split(".")[0]
    sam_file_path = os.path.join(sample_dir,f'{sample_name}.sam')
    bam_file_path = os.path.join(sample_dir,f'{sample_name}.bam')
    bam_tmp_path = bam_file_path + ".tmp"

    # Check if the output file already exists, if so skip the alignment.
    # Require a complete BAM (valid header + BGZF EOF block): a truncated BAM
    # from an interrupted `samtools view` must be re-run, not reused.
    if bam_ok(bam_file_path):
        log.info("Alignment has already been executed. Skip the alignment.")
        return bam_file_path

    # Remove any partial leftovers from a previous failed run so we start clean.
    for p in (sam_file_path, bam_file_path, bam_tmp_path):
        if os.path.exists(p):
            os.remove(p)
    try:
        # Run alignment with default bwa-mem behaviour so supplementary/SA
        # evidence is preserved. No read group is added (no downstream tool
        # needs it in TLA V1).
        bwa_command = ["bwa", "mem", "-t", "16", reference_path, processed_fastq_file_1_path, processed_fastq_file_2_path]
        # Open sam output file
        with open(sam_file_path,"w") as output_file:
            # Run BWA
            subprocess.run(bwa_command, stdout=output_file, check = True)

        # Convert sam file to bam file needed for later analysis. Write to a
        # temp path and atomically rename only after a quickcheck passes.
        samtools_command = ["samtools","view","-@","16","-S","-b",sam_file_path]
        # Open bam output file
        with open(bam_tmp_path, "w") as output_file:
            # Run samtools
            subprocess.run(samtools_command, stdout=output_file, check=True)

        if not commit(bam_tmp_path, bam_file_path, validator=bam_ok):
            raise subprocess.CalledProcessError(1, "samtools view")
        return bam_file_path

    except subprocess.CalledProcessError:
        log.error("Alignment failed.")
        for p in (sam_file_path, bam_file_path, bam_tmp_path):
            if os.path.exists(p):
                os.remove(p)
        return False


def _run_markdup_pipeline(input_bam, output_bam):
    """Stream ``input_bam`` through ``sort -n | fixmate -m | sort | markdup``.

    The chain is name-sort -> fixmate -> coordinate-sort -> markdup so
    duplicate flags (0x400) are set correctly. No intermediate BAMs hit the
    disk; only the final coordinate-sorted, duplicate-marked BAM is written to
    ``output_bam``. Raises ``CalledProcessError`` if any stage exits non-zero.
    """
    threads = "16"
    namesort = subprocess.Popen(
        ["samtools", "sort", "-n", "-@", threads, "-O", "bam", input_bam],
        stdout=subprocess.PIPE)
    fixmate = subprocess.Popen(
        ["samtools", "fixmate", "-m", "-@", threads, "-O", "bam", "-", "-"],
        stdin=namesort.stdout, stdout=subprocess.PIPE)
    namesort.stdout.close()  # allow namesort to receive SIGPIPE if fixmate exits
    coordsort = subprocess.Popen(
        ["samtools", "sort", "-@", threads, "-O", "bam", "-"],
        stdin=fixmate.stdout, stdout=subprocess.PIPE)
    fixmate.stdout.close()
    markdup = subprocess.Popen(
        ["samtools", "markdup", "-@", threads, "-", output_bam],
        stdin=coordsort.stdout)
    coordsort.stdout.close()

    markdup.communicate()  # block until the tail of the pipeline finishes
    stages = (("samtools sort -n", namesort), ("samtools fixmate", fixmate),
              ("samtools sort", coordsort), ("samtools markdup", markdup))
    for name, proc in stages:
        code = proc.wait()
        if code != 0:
            raise subprocess.CalledProcessError(code, name)


def samtools_runner(bam_file_path, flagstat_path):
    """Coordinate-sort + mark duplicates + index a BAM, write flagstat summary.

    Name-sorts the BAM, runs ``fixmate -m`` -> coordinate-sort -> ``markdup`` so
    PCR/optical duplicates are flagged (0x400) but kept. Writes a sorted BAM +
    ``.bai`` and a TSV flagstat summary (``flagstat_path``). Output dirs must
    already exist.

    Returns (sorted_bam, flagstat) on success, or (False, False) on failure.
    """
    sorted_bam_file_path = bam_file_path.replace(".bam", ".sorted.bam")
    indexed_bam_file_path = sorted_bam_file_path + ".bai"
    flagstat_file_path = flagstat_path

    # No-clobber: skip only if every output is present and the sorted BAM
    # passes a quickcheck (valid header + BGZF EOF block).
    if bam_ok(sorted_bam_file_path) and os.path.exists(indexed_bam_file_path) \
            and os.path.exists(flagstat_file_path):
        log.info("Samtools has already been executed on file {}".format(
            os.path.basename(bam_file_path)))
        return sorted_bam_file_path, flagstat_file_path

    sorted_tmp = sorted_bam_file_path + ".tmp"
    indexed_tmp = indexed_bam_file_path + ".tmp"
    flagstat_tmp = flagstat_file_path + ".tmp"

    # Remove any partial leftovers from a previous failed run so we start clean.
    for p in (sorted_bam_file_path, indexed_bam_file_path, flagstat_file_path,
              sorted_tmp, indexed_tmp, flagstat_tmp):
        if os.path.exists(p):
            os.remove(p)

    try:
        log.info("Coordinate-sorting and marking duplicates on the .BAM file produced.")
        _run_markdup_pipeline(bam_file_path, sorted_tmp)
        if not commit(sorted_tmp, sorted_bam_file_path, validator=bam_ok):
            raise subprocess.CalledProcessError(1, "samtools markdup")

        log.info("Indexing the sorted .BAM file.")
        subprocess.check_call(["samtools", "index", "-@", "16", sorted_bam_file_path, indexed_tmp])
        if not commit(indexed_tmp, indexed_bam_file_path, validator=is_nonempty):
            raise subprocess.CalledProcessError(1, "samtools index")

        log.info("Writing flagstat alignment summary (TSV).")
        with open(flagstat_tmp, "w") as fh:
            subprocess.check_call(
                ["samtools", "flagstat", "-O", "tsv", sorted_bam_file_path], stdout=fh)
        os.replace(flagstat_tmp, flagstat_file_path)

        return sorted_bam_file_path, flagstat_file_path

    except subprocess.CalledProcessError as e:
        log.error("SAMtools processing failed: {}".format(e))
        for p in (sorted_bam_file_path, indexed_bam_file_path, flagstat_file_path,
                  sorted_tmp, indexed_tmp, flagstat_tmp):
            if os.path.exists(p):
                os.remove(p)
        return False, False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Align paired FASTQ to a (mixed) reference and produce a "
                    "sorted+indexed BAM plus a flagstat summary.")
    parser.add_argument("-1", "--fastq1", required=True, help="R1 FASTQ")
    parser.add_argument("-2", "--fastq2", required=True, help="R2 FASTQ")
    parser.add_argument("-r", "--reference", required=True, help="Indexed (mixed) reference FASTA")
    parser.add_argument("-o", "--out_dir", required=True, help="Directory for BAM outputs")
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    bam = alignment(args.fastq1, args.fastq2, args.reference, args.out_dir)
    flagstat_path = bam.replace(".bam", ".flagstats.tsv")
    sorted_bam, flagstat = samtools_runner(bam, flagstat_path)
    log.info(f"sorted BAM: {sorted_bam}")
    log.info(f"flagstat:   {flagstat}")