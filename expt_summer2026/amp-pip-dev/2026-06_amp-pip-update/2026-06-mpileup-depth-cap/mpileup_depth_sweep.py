# 2026-06-mpileup-depth-cap/mpileup_depth_sweep.py
# Diagnostic: sweep samtools mpileup -d (max depth) and measure its effect on
# the depth reported in VarScan consensus output. Reuses existing sorted BAMs;
# does not re-align and does not touch production outputs.
# Usage: conda activate variant-calling; python mpileup_depth_sweep.py <amp-pip-dev>/runs/T3_rebuild
import glob
import os
import subprocess
import sys

# analysis.py lives in the sibling 2026-06-ref_coord_map experiment
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "2026-06-ref_coord_map"))

import plotly.graph_objects as go

from analysis import read_consensus, pad_to_reference

D_VALUES = [0, 8000, 12000, 16000]


def max_consensus_depth(consensus_path):
    """Max (Reads1 + Reads2) over all positions in a VarScan consensus file."""
    max_depth = 0
    with open(consensus_path) as fh:
        next(fh)  # header
        for line in fh:
            cols = line.split("\t")
            if len(cols) < 6 or not cols[1].isdigit():
                continue
            max_depth = max(max_depth, int(cols[4]) + int(cols[5]))
    return max_depth


def fasta_seq_name(fa_path):
    """First sequence name (first whitespace token after '>') in a FASTA file."""
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                return line[1:].split()[0]
    return None


def fasta_length(fa_path):
    """Total sequence length (bp) across all records in a FASTA file."""
    length = 0
    with open(fa_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                length += len(line.strip())
    return length


def plot_consensus_coverage(consensus_path, reference_length, sample, d, out_jpeg):
    """Plot WT vs variant coverage from one -d value's consensus file, padded to
    the full reference length, and save it as a JPEG."""
    positions, wt, variant = read_consensus(consensus_path)
    positions, (wt, variant) = pad_to_reference(positions, [wt, variant], reference_length)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=positions, y=wt, mode="lines", name="WT (reference)"))
    fig.add_trace(go.Scatter(x=positions, y=variant, mode="lines", name="Variant"))
    fig.update_layout(
        title="Consensus coverage for {} (mpileup -d {})".format(sample, d),
        xaxis_title="Position", yaxis_title="Number of reads",
    )
    fig.write_image(out_jpeg)


def bam_contigs(bam_path):
    """Set of @SQ SN: contig names in a BAM header."""
    header = subprocess.check_output(["samtools", "view", "-H", bam_path]).decode()
    contigs = set()
    for line in header.splitlines():
        if line.startswith("@SQ"):
            for field in line.split("\t"):
                if field.startswith("SN:"):
                    contigs.add(field[3:])
    return contigs


def select_reference(run_dir, bams):
    """Pick the reference whose sequence name matches a BAM contig.

    Falls back to the first glob match when no FASTA name matches (preserves the
    original behavior for runs with a single reference).
    """
    candidates = sorted(glob.glob(os.path.join(run_dir, "reference", "*.fa")))
    contigs = set()
    for bam in bams:
        contigs |= bam_contigs(bam)
    for fa in candidates:
        if fasta_seq_name(fa) in contigs:
            return fa
    return candidates[0]


def main():
    run_dir = sys.argv[1]
    bams = sorted(glob.glob(os.path.join(run_dir, "Results", "*", "*.sorted.bam")))
    reference = select_reference(run_dir, bams)
    reference_length = fasta_length(reference)
    print("reference: {}".format(reference))
    print("sample\td\tmax_consensus_depth")
    rows = []
    for bam in bams:
        sample_dir = os.path.dirname(bam)
        sample = os.path.basename(sample_dir)
        out_dir = os.path.join(sample_dir, "depth-sweep")
        os.makedirs(out_dir, exist_ok=True)
        for d in D_VALUES:
            mpileup_path = os.path.join(out_dir, "{}.d{}.mpileup".format(sample, d))
            consensus_path = os.path.join(out_dir, "{}.d{}.consensus.vcf".format(sample, d))
            subprocess.check_call(
                "samtools mpileup -f {} {} -A -d {} > {}".format(reference, bam, d, mpileup_path),
                shell=True,
            )
            subprocess.check_call(
                "varscan pileup2cns {} > {}".format(mpileup_path, consensus_path),
                shell=True,
            )
            md = max_consensus_depth(consensus_path)
            rows.append((sample, d, md))
            plot_path = os.path.join(out_dir, "{}.d{}.consensus_coverage.jpeg".format(sample, d))
            plot_consensus_coverage(consensus_path, reference_length, sample, d, plot_path)
            print("{}\t{}\t{}".format(sample, d, md))

    print("\n=== summary (max consensus depth) ===")
    print("d\\sample\t" + "\t".join(sorted({r[0] for r in rows})))
    samples = sorted({r[0] for r in rows})
    by = {(s, d): m for s, d, m in rows}
    for d in D_VALUES:
        print(str(d) + "\t" + "\t".join(str(by.get((s, d), "")) for s in samples))


if __name__ == "__main__":
    main()
