#!/usr/bin/env python3

"""Pure tool wrappers for plasmid/vector variant calling (TLA Stage 5, Part 1).

Active callers used by vector_variants.py:

    - ``plasmid_mpileup``  -- samtools mpileup restricted to the plasmid contig
    - ``varscan_call``     -- VarScan pileup2snp / pindel2indel / pileup2cns

Reserved (not yet wired into the pipeline; kept for future Pindel integration):

    - ``get_insert_size``  -- mean insert size from ``samtools stats`` (region)
    - ``filter_deletions`` -- filter pindel2vcf _D output to clean DEL records
    - ``pindel_call``      -- Pindel restricted to the plasmid contig

Adapted from the amp-pip-dev wrappers; rewritten here as argument-driven,
no-clobber helpers with no hard-coded paths or sample IDs.
"""

import os
import subprocess

# VarScan defaults tuned for a non-diploid plasmid: keep the coverage floor low
# so edge/low-depth positions are not dropped, and do the real filtering on
# variant allele frequency (a plasmid has variants at arbitrary fractions, not
# 0/1 / 1/1 genotypes).
DEFAULT_MIN_COVERAGE = 1500
DEFAULT_MIN_VAR_FREQ = 0.01

# Pindel sensitivity filter (fraction of reads supporting the event); matches
# the amp-pip-dev setting.
PINDEL_MIN_SUPPORT = "0.9"
PINDEL_THREADS = "16"
PINDEL2VCF_FILTER = "30"  # -e: minimum supporting reads for a called event


def _nonempty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


def plasmid_mpileup(sorted_bam, mixed_fa, plasmid_contig, out_path):
    """Write an mpileup for ``plasmid_contig`` only, uncapped depth.

    ``samtools mpileup -f <ref> -r <contig> -A -d 0`` -- ``-d 0`` removes the
    default 8000x depth cap (plasmid pile-ups in TLA are very deep) and ``-A``
    keeps anomalous read pairs. No-clobber + atomic (``.tmp`` renamed on
    success). Returns the mpileup path.
    """
    if _nonempty(out_path):
        return out_path
    tmp = out_path + ".tmp"
    if os.path.exists(tmp):
        os.remove(tmp)
    try:
        with open(tmp, "w") as fh:
            subprocess.check_call(
                ["samtools", "mpileup", "-f", mixed_fa,
                 "-r", plasmid_contig, "-A", "-d", "0", sorted_bam],
                stdout=fh)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    os.replace(tmp, out_path)
    return out_path


def varscan_call(mpileup_path, out_dir, prefix,
                 min_coverage=DEFAULT_MIN_COVERAGE,
                 min_var_freq=DEFAULT_MIN_VAR_FREQ):
    """Run VarScan SNP / indel / consensus on an mpileup.

    Returns ``(snp_path, indel_path, consensus_path)``. Frequency-aware:
    ``--min-coverage`` is a low floor and ``--min-var-freq`` is the real filter.
    No-clobber: if all three outputs already exist they are returned as-is.
    Output dir must already exist (created by the caller).
    """
    snp = os.path.join(out_dir, prefix + ".varscan.snp.tsv")
    indel = os.path.join(out_dir, prefix + ".varscan.indel.tsv")
    cns = os.path.join(out_dir, prefix + ".varscan.cns.tsv")

    if _nonempty(snp) and _nonempty(indel) and _nonempty(cns):
        return snp, indel, cns

    common = ["--min-coverage", str(min_coverage),
              "--min-var-freq", str(min_var_freq)]
    for sub, out_path in (("pileup2snp", snp),
                          ("pileup2indel", indel),
                          ("pileup2cns", cns)):
        tmp = out_path + ".tmp"
        if os.path.exists(tmp):
            os.remove(tmp)
        try:
            with open(tmp, "w") as fh:
                subprocess.check_call(
                    ["varscan", sub, mpileup_path] + common, stdout=fh)
        except BaseException:
            if os.path.exists(tmp):
                os.remove(tmp)
            raise
        os.replace(tmp, out_path)
    return snp, indel, cns


# ---------------------------------------------------------------------------
# Reserved functions (not wired into the pipeline; kept for future Pindel work)
# ---------------------------------------------------------------------------

def get_insert_size(sorted_bam, plasmid_contig):
    """Return the mean insert size (rounded int) from ``samtools stats``.

    Restricted to ``plasmid_contig`` so the value reflects plasmid-mapped pairs,
    not the whole (mostly host) library. Reads the ``SN  insert size average:``
    line. Raises ``ValueError`` if it is absent or non-positive (too few
    properly-paired plasmid reads to estimate).
    """
    out = subprocess.check_output(
        ["samtools", "stats", sorted_bam, plasmid_contig], text=True)
    for line in out.splitlines():
        if line.startswith("SN") and "insert size average:" in line:
            value = float(line.split("insert size average:")[1].strip().split()[0])
            if value <= 0:
                raise ValueError(
                    f"Non-positive mean insert size ({value}) for "
                    f"{plasmid_contig}; too few properly-paired plasmid reads?")
            return int(round(value))
    raise ValueError(
        f"'insert size average' not found in samtools stats for {plasmid_contig}")


def filter_deletions(deletion_vcf_path):
    """Keep only clean ``DEL`` records (``HOMLEN=0``) from a pindel2vcf ``_D`` VCF.

    Reserved -- not called by the active pipeline.
    """
    if not _nonempty(deletion_vcf_path):
        return None
    with open(deletion_vcf_path) as fh:
        rows = [ln.split("\t") for ln in fh.read().split("\n")
                if ln and not ln.startswith("#")]
    kept = [r for r in rows
            if len(r) > 7 and "DEL" in r[7] and "HOMLEN=0" in r[7]]
    filtered_path = deletion_vcf_path.replace("_D", "_D_Filtered")
    with open(filtered_path, "w") as fh:
        for r in kept:
            fh.write("\t".join(r) + "\n")
    return filtered_path


def _pindel_intermittent_files(pindel_out_prefix):
    """Return the SI/LI/D Pindel output files that exist for this prefix."""
    suffixes = ["_SI", "_LI", "_D"]
    return [pindel_out_prefix + s for s in suffixes
            if os.path.exists(pindel_out_prefix + s)]


def pindel_call(sorted_bam, mixed_fa, out_dir, prefix, insert_size,
                plasmid_contig):
    """Run Pindel restricted to the plasmid contig, then pindel2vcf.

    Optional caller (the stage leaves this off by default; DELLY covers the
    mid-size indel/SV range). Restricts to ``plasmid_contig`` with ``-c`` and
    uses ``insert_size`` (derived from ``get_insert_size``) in the config.
    Produces per-type VCFs (SI/LI/D) and a homology-filtered deletion VCF.
    Returns the list of VCF paths written. Output dir must already exist.
    """
    config_path = os.path.join(out_dir, prefix + ".pindel.config.txt")
    with open(config_path, "w") as fh:
        fh.write(f"{sorted_bam}\t{insert_size}\t{prefix}\n")

    pindel_out_prefix = os.path.join(out_dir, prefix + ".pindel")
    if not os.path.exists(pindel_out_prefix + "_SI"):
        subprocess.check_call(
            ["pindel", "-f", mixed_fa, "-i", config_path,
             "-c", plasmid_contig, "-E", PINDEL_MIN_SUPPORT,
             "-T", PINDEL_THREADS, "-o", pindel_out_prefix])

    vcf_paths = []
    for intermittent in _pindel_intermittent_files(pindel_out_prefix):
        vcf_path = intermittent + f".{PINDEL2VCF_FILTER}.vcf"
        if not _nonempty(vcf_path):
            subprocess.check_call(
                ["pindel2vcf", "-p", intermittent, "-r", mixed_fa,
                 "-R", "x", "-d", "00000000",
                 "-e", PINDEL2VCF_FILTER, "-v", vcf_path])
        vcf_paths.append(vcf_path)
        if intermittent.endswith("_D"):
            filtered = filter_deletions(vcf_path)
            if filtered:
                vcf_paths.append(filtered)
    return vcf_paths


if __name__ == "__main__":
    pass
