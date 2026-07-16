#!/usr/bin/env python
"""GATK HaplotypeCaller wrapper + per-variant read-count table.

All GATK logic lives here:
  - gatk_runner()  : Python port of the former gatk_variant_calling.sh.
  - run_gatk()     : orchestrator used by the pipeline (calls gatk_runner,
                     then writes the per-variant read-count table).
  - count helpers  : parse the VCF + BAM into gatk_variant.tsv.

Original bash script created by Lilin Huang.
Updated by Cindy Ren 2026-06-02 (in the bash version, carried over here):
  - No Picard AddOrReplaceReadGroups / re-sort: the input BAM is already
    coordinate-sorted and carries @RG from BWA alignment (see align.py).
  - Skip CreateSequenceDictionary when the .dict already exists, since the
    reference directory may be a read-only mount.
  - --max-reads-per-alignment-start = 2 * mapped reads (disable downsampling).

Requires: gatk (>=4.4), samtools, pysam.
"""

import os
import subprocess

import pysam

from pipeline_io import gzip_ok


def gatk_runner(reference_path, sorted_bam_path, output_dir):
    """Run GATK HaplotypeCaller; return the bgzipped VCF path (or False on failure)."""
    if not os.access(reference_path, os.R_OK):
        print(f"[Takeda Pipeline Message] Reference not readable: {reference_path}")
        return False
    if not os.access(sorted_bam_path, os.R_OK):
        print(f"[Takeda Pipeline Message] BAM not readable: {sorted_bam_path}")
        return False

    os.makedirs(output_dir, exist_ok=True)
    output_vcf_path = os.path.join(output_dir, "gatk_variants.vcf.gz")
    tmp_vcf_path = os.path.join(output_dir, "gatk_variants.tmp.vcf.gz")

    # Require a complete bgzipped VCF (not just present) so a truncated file
    # from an interrupted/killed HaplotypeCaller run is re-run, not reused.
    if gzip_ok(output_vcf_path):
        print("[Takeda Pipeline Message] GATK HaplotypeCaller already executed; skipping.")
        return output_vcf_path

    # Remove any partial leftovers (final + temp VCF + tabix index) from a
    # failed run so we start clean.
    for p in (output_vcf_path, output_vcf_path + ".tbi",
              tmp_vcf_path, tmp_vcf_path + ".tbi"):
        if os.path.exists(p):
            os.remove(p)

    try:
        # Number of mapped reads (used to size --max-reads-per-alignment-start).
        number_of_reads = int(subprocess.check_output(
            ["samtools", "view", "-c", "-F", "4", sorted_bam_path]).strip())
        print(f"[Takeda Pipeline Message] Number of mapped reads: {number_of_reads}")

        # Create reference .dict next to the FASTA if missing. GATK accepts
        # either <ref>.dict (extension replaced) or <ref>.fa.dict.
        ref_no_ext = os.path.splitext(reference_path)[0]
        if not os.path.exists(ref_no_ext + ".dict") and not os.path.exists(reference_path + ".dict"):
            subprocess.check_call(["gatk", "CreateSequenceDictionary", "-R", reference_path])

        # Write to a temp VCF and atomically rename only after a gzip integrity
        # check passes, so a killed HaplotypeCaller never leaves a partial VCF
        # at the canonical path.
        reads_consider = number_of_reads * 2
        subprocess.check_call([
            "gatk", "HaplotypeCaller",
            "-R", reference_path,
            "-I", sorted_bam_path,
            "--max-reads-per-alignment-start", str(reads_consider),
            "--stand-call-conf", "20",
            "--active-probability-threshold", "0.0002",
            "--min-base-quality-score", "5",
            "-O", tmp_vcf_path,
        ])
        if not gzip_ok(tmp_vcf_path):
            raise subprocess.CalledProcessError(1, "gatk HaplotypeCaller")
        os.replace(tmp_vcf_path, output_vcf_path)
        if os.path.exists(tmp_vcf_path + ".tbi"):
            os.replace(tmp_vcf_path + ".tbi", output_vcf_path + ".tbi")
    except subprocess.CalledProcessError as error:
        print(f"[Takeda Pipeline Message] GATK HaplotypeCaller failed: {error}")
        for p in (output_vcf_path, output_vcf_path + ".tbi",
                  tmp_vcf_path, tmp_vcf_path + ".tbi"):
            if os.path.exists(p):
                os.remove(p)
        return False

    print("[Takeda Pipeline Message] HaplotypeCaller executed successfully.")
    print(f"[Takeda Pipeline Message] Check the results in {output_dir}.")
    return output_vcf_path


def run_gatk(reference_path, sorted_bam_path, output_dir):
    """Run HaplotypeCaller, then write the per-variant read-count table."""
    output_dir = str(output_dir)
    gatk_vcf = gatk_runner(reference_path, sorted_bam_path, output_dir)
    if not gatk_vcf:
        return False
    counts_path = os.path.join(output_dir, "gatk_variant.tsv")
    if not os.path.exists(counts_path):
        count_reads_gatk(sorted_bam_path, gatk_vcf, counts_path)
    return gatk_vcf


def get_mapped_reads_count(bam_file_path):
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")
    mapped_reads = bamfile.mapped
    bamfile.close()
    return mapped_reads


def get_gatk_counts(vcf_file_path):
    counts = []
    with pysam.VariantFile(vcf_file_path) as vcf:
        for record in vcf:
            svtype = "DEL" if len(record.ref) > len(record.alts[0]) else "INS"

            start_pos = record.pos
            end_pos = record.pos + len(record.ref) - 1
            ad_values = record.samples[0]['AD']
            reads_count = ad_values[1] if len(ad_values) > 1 else 0
            wt_reads_count = ad_values[0] if len(ad_values) > 1 else 0
            sum_ad = sum(ad_values)
            dp_format = record.samples[0].get('DP', None)
            dp_info = record.info.get('DP', None)
            # pct_alt: ALT fraction among GATK-attributed reads (reads_count / sum_AD).
            pct_alt = round(100.0 * reads_count / sum_ad, 2) if sum_ad else None
            # pct_sum_ad: how much of the raw locus depth GATK could attribute to an allele (sum_AD / DP_INFO).
            pct_sum_ad = round(100.0 * sum_ad / dp_info, 2) if dp_info else None
            counts.append((start_pos + 1, end_pos, reads_count, svtype, wt_reads_count,
                           sum_ad, dp_format, dp_info, pct_alt, pct_sum_ad))
    return counts


def write_counts_to_file(output_file_path, total_mapped_reads, gatk_counts):
    total_ko_reads = sum(count[2] for count in gatk_counts)
    with open(output_file_path, 'w') as f:
        f.write(f"# total number of analyzed reads: {total_mapped_reads}\n")
        f.write(f"# number of unanalyzed/uninformative reads: {total_mapped_reads - total_ko_reads}\n")
        f.write("start_pos\tend_pos\treads_count\tSVTYPE\tWT_reads_count\tsum_AD\tDP_FORMAT\tDP_INFO\tpct_alt\tpct_sum_ad\n")
        for count in gatk_counts:
            f.write(f"{count[0]}\t{count[1]}\t{count[2]}\t{count[3]}\t{count[4]}\t"
                    f"{count[5]}\t{count[6]}\t{count[7]}\t{count[8]}\t{count[9]}\n")


def count_reads_gatk(bam_file_path, gatk_vcf_file_path, output_file_path):
    total_mapped_reads = get_mapped_reads_count(bam_file_path)
    gatk_counts = get_gatk_counts(gatk_vcf_file_path)
    write_counts_to_file(output_file_path, total_mapped_reads, gatk_counts)


if __name__ == "__main__":
    pass
