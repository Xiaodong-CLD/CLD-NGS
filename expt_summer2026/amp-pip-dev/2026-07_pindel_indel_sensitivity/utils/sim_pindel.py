#!/usr/bin/env python3
"""Sim-specific Pindel runner: raw calls only, no production filtering.

A trimmed copy of the production pindel_wrapper.py tailored to the benchmark.
Deliberately different from the root wrapper (which must stay untouched):
  - runs pindel2vcf ONCE at -e 1 per SV type and stops there (the base VCF keeps
    every call with the AD supporting-read field intact, so evaluate.py applies
    every -e / fraction threshold post-hoc in Python -- one Pindel run, infinite
    re-scoring);
  - does NOT run the deletion HOMLEN/RPL filter (those microhomology-shifted
    calls are exactly the near-misses the benchmark wants to keep and inspect);
  - does NOT write the absolute/fractional ladder files or Summary.100.vcf.

Emits base per-type VCFs <prefix>_SI.vcf / <prefix>_LI.vcf / <prefix>_D.vcf in
the sorted BAM's pindel-output/ directory.
"""

import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config  # noqa: E402

# Pindel insert size is no longer hardcoded here; it is passed in as a
# parameter derived from config.amplicon_insert_size(amplicon_len).
# Type suffixes Pindel writes that we score (short/long insertions, deletions).
TYPE_SUFFIXES = ["_SI", "_LI", "_D"]


def _existing_type_files(prefix):
    return [prefix + s for s in TYPE_SUFFIXES if os.path.exists(prefix + s)]


def _chrom_name(ref_path):
    """Return the first FASTA record name from ref_path."""
    with open(ref_path) as fh:
        for line in fh:
            if line.startswith(">"):
                return line[1:].split()[0].strip()
    raise ValueError(f"No FASTA records in {ref_path}")


def _chrom_name_from_bam(bam_path):
    """Return the first reference sequence name from the BAM @SQ header.

    Used instead of _chrom_name(ref) when the BAM was aligned to a WT
    reference whose header name differs from the mutant reference FASTA
    (e.g. amplicon_LAGA_LC_WT vs amplicon_LAGA_LC_del_L10_p355).
    """
    result = subprocess.run(
        ["samtools", "view", "-H", bam_path],
        capture_output=True, text=True, check=True)
    for line in result.stdout.splitlines():
        if line.startswith("@SQ"):
            for field in line.split("\t"):
                if field.startswith("SN:"):
                    return field[3:]
    raise ValueError(f"No @SQ lines found in BAM header: {bam_path}")


def run_sim_pindel(sorted_bam_file_path, reference_file_path, insert_size):
    """Run Pindel + pindel2vcf (-e 1) for one condition.

    insert_size: Pindel's -i field (fragment length between adapters); use
    config.amplicon_insert_size(amplicon_len) to derive from the dataset.

    Returns the list of base per-type VCF paths written.
    """
    bam_dir = os.path.dirname(sorted_bam_file_path)
    bam_name = os.path.basename(sorted_bam_file_path)
    out_dir = os.path.join(bam_dir, "pindel-output")
    os.makedirs(out_dir, exist_ok=True)

    config_file_path = os.path.join(
        out_dir, bam_name.replace(".bam", ".pindel.config.txt"))
    prefix = os.path.join(out_dir, bam_name.replace(".bam", ".pindel"))

    with open(config_file_path, "w") as fh:
        fh.write(f"{sorted_bam_file_path}\t{insert_size}\t{bam_name}")

    # Run Pindel to produce the intermittent per-type files (skip if present).
    if not os.path.exists(prefix + "_SI"):
        # Use the chromosome name from the BAM header, not the mutant reference
        # FASTA. The BAM is always aligned to the WT reference, so its @SQ SN
        # is the authoritative name (e.g. amplicon_LAGA_LC_WT), which may differ
        # from the mutant FASTA header (e.g. amplicon_LAGA_LC_del_L10_p355).
        chrom = _chrom_name_from_bam(sorted_bam_file_path)
        cmd = (f"pindel -f {reference_file_path} -i {config_file_path} "
               f"-E 0.9 -T 1 -c {chrom} -o {prefix}")
        print(f"[sim_pindel] {cmd}")
        if subprocess.call(cmd, shell=True) != 0:
            raise RuntimeError(f"pindel failed on {bam_name}")
    else:
        print(f"[sim_pindel] reusing existing Pindel intermittent files for "
              f"{bam_name}")

    # One pindel2vcf per type at -e 1 (keep every call; AD field carries support).
    base_vcfs = []
    for type_file in _existing_type_files(prefix):
        base_vcf = os.path.join(out_dir, os.path.basename(type_file) + ".vcf")
        if not os.path.exists(base_vcf):
            cmd = (f"pindel2vcf -p {type_file} -r {reference_file_path} "
                   f"-R x -d 00000000 -e {config.CEILING_E} -v {base_vcf}")
            print(f"[sim_pindel] {cmd}")
            if subprocess.call(cmd, shell=True) != 0:
                raise RuntimeError(f"pindel2vcf failed on {type_file}")
        base_vcfs.append(base_vcf)
    return base_vcfs
