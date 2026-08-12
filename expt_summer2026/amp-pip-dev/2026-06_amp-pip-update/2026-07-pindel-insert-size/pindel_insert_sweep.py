# 2026-07-pindel-insert-size/pindel_insert_sweep.py
# Diagnostic: sweep the Pindel config insert-size field (the hard-coded 398 in
# pindel_wrapper.py) and measure its effect on the indel calls that reach the
# *.100.vcf summary (SI + LI + filtered D). Also reports the empirical insert
# size and read length from each BAM so 398 can be judged against the data.
#
# Reuses existing sorted BAMs; does NOT re-align and does NOT touch production
# pindel-output/ (writes only under Results/<sample>/pindel-insert-sweep/).
#
# Usage: conda activate variant-calling; python pindel_insert_sweep.py <amp-pip-dev>/runs/<RUN_DIR>
import glob
import os
import subprocess
import sys

# Import the production D-filter so the experiment counts deletions the same way
# the pipeline summary does (SVTYPE=DEL & HOMLEN=0).
# Path: this file is at expt_summer2026/amp-pip-dev/2026-06_amp-pip-update/2026-07-pindel-insert-size/
# amp-pip-dev/ is 4 levels up.
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "amp-pip-dev"))
sys.path.insert(0, REPO_ROOT)
from filter_deletions import filter_deletions_from_pindel2vcf  # noqa: E402

# Fixed reference points bracketing the current literal (398). The per-sample
# empirical mean insert size is merged into this set at runtime.
BASE_INSERT_VALUES = [200, 300, 398, 500, 600]

# Production Pindel / pindel2vcf flags (keep in sync with pindel_wrapper.py).
PINDEL_FILTER = "-E 0.9"
PINDEL_THREADS = 16
PINDEL2VCF_E = 100
SUFFIXES_OF_INTEREST = ["_SI", "_LI", "_D"]


def samtools_stats_metrics(bam_path):
    """Return (read_len_avg, read_len_max, insert_avg, insert_sd) from samtools stats.

    Values Pindel cares about live in the SN (summary numbers) block:
      'average length', 'maximum length', 'insert size average',
      'insert size standard deviation'.
    """
    out = subprocess.check_output(["samtools", "stats", bam_path]).decode()
    metrics = {}
    for line in out.splitlines():
        if not line.startswith("SN\t"):
            continue
        # SN\t<key>:\t<value>\t...
        parts = line.split("\t")
        key = parts[1].rstrip(":")
        metrics[key] = parts[2]
    return (
        float(metrics.get("average length", 0.0)),
        int(float(metrics.get("maximum length", 0.0))),
        float(metrics.get("insert size average", 0.0)),
        float(metrics.get("insert size standard deviation", 0.0)),
    )


def count_vcf_records(vcf_path):
    """Number of non-header (non-'#') records in a VCF file (0 if missing)."""
    if not os.path.exists(vcf_path):
        return 0
    n = 0
    with open(vcf_path) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                n += 1
    return n


def run_pindel_for_insert(bam_path, reference, sample, insert, out_dir):
    """Run Pindel + pindel2vcf for one insert value; return per-type call counts.

    Returns dict {'SI': n, 'LI': n, 'D': n, 'total': n}. Skips the Pindel step if
    its intermittent files already exist (resume-friendly, like production)."""
    os.makedirs(out_dir, exist_ok=True)
    config_path = os.path.join(out_dir, "{}.pindel.config.txt".format(sample))
    prefix = os.path.join(out_dir, "{}.pindel".format(sample))

    with open(config_path, "w") as fh:
        fh.write("{}\t{}\t{}".format(bam_path, insert, sample))

    if not os.path.exists(prefix + "_SI"):
        subprocess.call(
            "pindel -f {} -i {} {} -T {} -o {}".format(
                reference, config_path, PINDEL_FILTER, PINDEL_THREADS, prefix
            ),
            shell=True,
        )

    counts = {"SI": 0, "LI": 0, "D": 0}
    for suffix in SUFFIXES_OF_INTEREST:
        intermittent = prefix + suffix
        if not os.path.exists(intermittent):
            continue
        vcf_path = intermittent + ".{}.vcf".format(PINDEL2VCF_E)
        subprocess.call(
            "pindel2vcf -p {} -r {} -R x -d 00000000 -e {} -v {}".format(
                intermittent, reference, PINDEL2VCF_E, vcf_path
            ),
            shell=True,
        )
        if suffix == "_D":
            # Apply production DEL/HOMLEN=0 filter, then count the filtered file.
            filter_deletions_from_pindel2vcf(vcf_path)
            counts["D"] = count_vcf_records(vcf_path.replace("_D", "_D_Filtered"))
        elif suffix == "_SI":
            counts["SI"] = count_vcf_records(vcf_path)
        elif suffix == "_LI":
            counts["LI"] = count_vcf_records(vcf_path)

    counts["total"] = counts["SI"] + counts["LI"] + counts["D"]
    return counts


def select_reference(run_dir):
    """First reference FASTA in the run's reference/ directory."""
    candidates = sorted(glob.glob(os.path.join(run_dir, "reference", "*.fa")))
    if not candidates:
        raise SystemExit("No reference FASTA under {}/reference/".format(run_dir))
    return candidates[0]


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python pindel_insert_sweep.py <RUN_DIR>")
    run_dir = sys.argv[1]
    reference = select_reference(run_dir)
    bams = sorted(glob.glob(os.path.join(run_dir, "Results", "*", "*.sorted.bam")))
    if not bams:
        raise SystemExit("No *.sorted.bam under {}/Results/*/".format(run_dir))
    print("reference: {}".format(reference))

    for bam in bams:
        sample_dir = os.path.dirname(bam)
        sample = os.path.basename(sample_dir)
        rl_avg, rl_max, ins_avg, ins_sd = samtools_stats_metrics(bam)
        print("\n=== {} ===".format(sample))
        print("read length: avg={:.1f} max={}".format(rl_avg, rl_max))
        print("empirical insert size: mean={:.1f} sd={:.1f}".format(ins_avg, ins_sd))

        insert_values = sorted(set(BASE_INSERT_VALUES) | {round(ins_avg)})
        sweep_root = os.path.join(sample_dir, "pindel-insert-sweep")

        print("insert\tSI\tLI\tD_filt\ttotal")
        for insert in insert_values:
            out_dir = os.path.join(sweep_root, "insert-{}".format(insert))
            counts = run_pindel_for_insert(bam, reference, sample, insert, out_dir)
            tag = "  <- empirical" if insert == round(ins_avg) else (
                "  <- current" if insert == 398 else "")
            print("{}\t{}\t{}\t{}\t{}{}".format(
                insert, counts["SI"], counts["LI"], counts["D"], counts["total"], tag))


if __name__ == "__main__":
    main()
