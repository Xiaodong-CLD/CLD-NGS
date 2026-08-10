# Follow-up exploratory analysis for completed amp-pip TLA runs.
# Usage: python analysis.py <amp-pip-dev>/runs/P3_TLA
#
# For each sample in <run>/Results/<sample>/:
#   - replot coverage from <sample>.depth.txt with HC/LC/GS CDS regions shaded
#   - write <sample>.cds_depth.txt with per-feature depth blocks
# Assumes a completed, well-formed run folder (no I/O checks).

import glob
import os
import sys

from Bio import SeqIO
import plotly.graph_objects as go

# label substring -> (group key, color)
CDS_GROUPS = [
    ("_HC", "HC", "rgba(31,119,180,0.20)"),
    ("_LC", "LC", "rgba(44,160,44,0.20)"),
    ("Glutamine Synthetase", "GS", "rgba(214,39,40,0.20)"),
]


def get_cds_regions(gb_file):
    """Return [(group, name, start, end)] for the HC/LC/GS CDS, longest per group."""
    record = SeqIO.read(gb_file, "genbank")
    picked = {}
    for f in record.features:
        if f.type != "CDS":
            continue
        labels = f.qualifiers.get("label", []) + f.qualifiers.get("gene", [])
        start = int(f.location.start) + 1
        end = int(f.location.end)
        for needle, group, _color in CDS_GROUPS:
            if any(needle in lab for lab in labels):
                length = end - start
                if group not in picked or length > picked[group][3] - picked[group][2]:
                    picked[group] = (group, labels[0], start, end)
                break
    return [picked[g] for _, g, _ in CDS_GROUPS if g in picked]


def get_reference_length(gb_file):
    """Total length (bp) of the GenBank reference, used to pad coverage so plots
    span the full amplicon."""
    return len(SeqIO.read(gb_file, "genbank").seq)


def read_depth(depth_file):
    """Return (positions, depths) from a samtools depth file."""
    positions, depths = [], []
    with open(depth_file) as fh:
        for line in fh:
            _contig, pos, depth = line.split("\t")
            positions.append(int(pos))
            depths.append(int(depth))
    return positions, depths


def read_consensus(consensus_file):
    """Return (positions, wt, variant) read counts from a VarScan consensus file.

    wt      = Reads1Plus + Reads1Minus (cols 14+15): reads matching reference.
    variant = Reads2Plus + Reads2Minus (cols 16+17): reads supporting the variant.
    These are quality-filtered counts (min-coverage 1500), not raw depth.
    """
    positions, wt, variant = [], [], []
    with open(consensus_file) as fh:
        next(fh)  # header row
        for line in fh:
            if not line.strip():
                continue
            cols = line.split("\t")
            positions.append(int(cols[1]))
            wt.append(int(cols[14]) + int(cols[15]))
            variant.append(int(cols[16]) + int(cols[17]))
    return positions, wt, variant


def pad_to_reference(positions, series, reference_length):
    """Pad coverage with zero-coverage flanks so the plot spans the full reference:
    positions before the first covered base and after the last, up to
    reference_length. ``series`` is a list of value-lists aligned with positions."""
    if not positions:
        return positions, series
    order = sorted(range(len(positions)), key=lambda i: positions[i])
    positions = [positions[i] for i in order]
    series = [[s[i] for i in order] for s in series]
    lead = list(range(1, positions[0]))
    tail = list(range(positions[-1] + 1, reference_length + 1))
    positions = lead + positions + tail
    series = [[0] * len(lead) + s + [0] * len(tail) for s in series]
    return positions, series


def add_cds_shading(fig, regions):
    """Shade HC/LC/GS CDS regions on a figure using the group colors."""
    color_by_group = {g: c for _, g, c in CDS_GROUPS}
    for group, name, start, end in regions:
        fig.add_vrect(
            x0=start, x1=end,
            fillcolor=color_by_group[group], line_width=0, layer="below",
            annotation_text=group, annotation_position="top left",
        )


# plot the annotated coverage using mpileup depth
def plot_annotated_coverage(positions, depths, regions, sample, out_jpeg):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=positions, y=depths, mode="lines", name="Depth"))
    add_cds_shading(fig, regions)
    fig.update_layout(
        title="Annotated coverage for " + sample,
        xaxis_title="Position", yaxis_title="Depth",
    )
    fig.write_image(out_jpeg)

# plot annotated WT/variant coverage using varscan consensus
def plot_wt_variant(positions, wt, variant, regions, sample, out_jpeg):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=positions, y=wt, mode="lines", name="WT (reference)"))
    fig.add_trace(go.Scatter(x=positions, y=variant, mode="lines", name="Variant"))
    add_cds_shading(fig, regions)
    fig.update_layout(
        title="WT vs variant coverage for " + sample,
        xaxis_title="Position", yaxis_title="Number of reads",
    )
    fig.write_image(out_jpeg)


def write_cds_depth(positions, depths, regions, out_txt):
    depth_by_pos = dict(zip(positions, depths))
    with open(out_txt, "w") as out:
        for group, name, start, end in regions:
            out.write("Feature: {} ({}); position: {}, {}\n".format(name, group, start, end))
            for pos in range(start, end + 1):
                out.write("{}\t{}\n".format(pos, depth_by_pos[pos]))
            out.write("\n")


def main():
    run_dir = sys.argv[1]
    gb_file = glob.glob(os.path.join(run_dir, "reference", "*.gb"))[0]
    regions = get_cds_regions(gb_file)
    reference_length = get_reference_length(gb_file)

    for depth_file in glob.glob(os.path.join(run_dir, "Results", "*", "*.depth.txt")):
        sample_dir = os.path.dirname(depth_file)
        sample = os.path.basename(sample_dir)
        positions, depths = read_depth(depth_file)

        plot_dir = os.path.join(sample_dir, "coverage-plot-output")
        os.makedirs(plot_dir, exist_ok=True)
        out_jpeg = os.path.join(plot_dir, sample + ".annotated_coverage.jpeg")
        plot_pos, (plot_depths,) = pad_to_reference(positions, [depths], reference_length)
        plot_annotated_coverage(plot_pos, plot_depths, regions, sample, out_jpeg)

        consensus_file = os.path.join(sample_dir, "varscan-output", sample + ".consensus.vcf")
        cons_pos, wt, variant = read_consensus(consensus_file)
        cons_pos, (wt, variant) = pad_to_reference(cons_pos, [wt, variant], reference_length)
        out_var_jpeg = os.path.join(plot_dir, sample + ".wt_variant_coverage.jpeg")
        plot_wt_variant(cons_pos, wt, variant, regions, sample, out_var_jpeg)

        out_txt = os.path.join(sample_dir, sample + ".cds_depth.txt")
        write_cds_depth(positions, depths, regions, out_txt)

        print("[analysis] {}: wrote {}, {} and {}".format(sample, out_jpeg, out_var_jpeg, out_txt))


if __name__ == "__main__":
    main()