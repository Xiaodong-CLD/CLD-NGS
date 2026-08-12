#!/usr/bin/env python3

"""TLA-seq -- per-base vector (plasmid) coverage plot (standalone).

From a sorted BAM (mixed host+vector reference) + the plasmid contig name:

    -> per-base depth over the plasmid   (<prefix>.vector_depth.tsv)
    -> linear per-base coverage plot      (<prefix>.vector_coverage.png)

Unlike the host coverage plot (100 kb-binned, multi-chromosome lanes), the
plasmid is small (~10 kb) so this draws depth at single-base resolution. The
x-axis always spans the full plasmid reference -- base 1 to the contig length
from the ``.fai`` -- so the plot starts at the first base and ends at the last
base regardless of where read coverage begins or ends. The y-axis is linear.

Runs standalone (CLI below) and as a helper of the Part 1 variant stage.
"""

import argparse
import os
import subprocess

import numpy as np
import plotly.graph_objects as go

from pipeline_io import is_nonempty, bam_ok, get_logger

log = get_logger()


def get_contig_length(fai_path, contig):
    """Return the length (bp) of ``contig`` from a ``.fai`` index."""
    with open(fai_path) as fh:
        for line in fh:
            fields = line.split("\t")
            if fields[0] == contig:
                return int(fields[1])
    raise ValueError(f"Contig {contig!r} not found in {fai_path}")


def read_depth(depth_path, length):
    """Read ``samtools depth`` output into a per-base array of size ``length``.

    ``samtools depth`` rows are ``contig<TAB>pos<TAB>depth`` (pos is 1-based).
    Positions absent from the file are left at 0. Returns a ``float`` array
    indexed 0..length-1 (position p -> index p-1).
    """
    depth = np.zeros(length, dtype=float)
    with open(depth_path) as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            pos = int(parts[1])
            if 1 <= pos <= length:
                depth[pos - 1] = float(parts[2])
    return depth


def plasmid_depth(sorted_bam, plasmid_contig, out_path):
    """Write per-base depth (all positions) for ``plasmid_contig`` to ``out_path``.

    ``samtools depth -a -r <contig>`` emits every base (``-a`` includes
    zero-coverage positions). No-clobber + atomic. Returns ``out_path``.
    """
    if is_nonempty(out_path):
        log.info(f"Vector depth already present, skipping: {out_path}")
        return out_path
    if not bam_ok(sorted_bam):
        raise FileNotFoundError(f"Sorted BAM missing or truncated: {sorted_bam}")
    tmp = out_path + ".tmp"
    if os.path.exists(tmp):
        os.remove(tmp)
    try:
        with open(tmp, "w") as fh:
            subprocess.check_call(
                ["samtools", "depth", "-a", "-r", plasmid_contig, sorted_bam],
                stdout=fh)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    os.replace(tmp, out_path)
    return out_path


def plot_vector_coverage(depth, length, out_png, label, contig, gridss_rows=None):
    """Draw a linear per-base coverage plot spanning the full plasmid reference.

    ``depth`` is a per-base array (index 0 -> position 1). The x-axis is fixed
    to ``[1, length]``. When ``gridss_rows`` is provided, vector–vector GRIDSS
    breakpoint endpoints on this contig are overlaid as vertical markers colored
    by NlaIII classification with a legend. Exported as PNG via kaleido 0.x
    (bundled Chromium; no system Chrome needed). Returns ``(1, length)``.
    """
    _NLAIII_COLORS = {
        "both_ends_near_nlaiii":   "crimson",
        "one_end_near_nlaiii":     "darkorange",
        "neither_end_near_nlaiii": "dimgray",
    }
    _DEFAULT_BP_COLOR = "darkviolet"

    positions = list(range(1, length + 1))
    y_max = float(np.max(depth)) * 1.1 if depth.size and np.max(depth) > 0 else 10
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=positions, y=list(depth),
        mode="lines",
        line=dict(color="steelblue", width=1.5),
        name="depth",
        showlegend=False,
    ))

    if gridss_rows:
        seen = {}  # category label -> color
        shapes = []
        for row in gridss_rows:
            cat = row.get("nlaiii_category", "")
            color = _NLAIII_COLORS.get(cat, _DEFAULT_BP_COLOR)
            cat_lbl = cat if cat else "vector_sv"
            seen.setdefault(cat_lbl, color)
            for pos_key, ctg_key in (("pos_a", "contig_a"), ("pos_b", "contig_b")):
                if row.get(ctg_key) == contig:
                    try:
                        p = int(row[pos_key])
                    except (ValueError, TypeError):
                        continue
                    if 1 <= p <= length:
                        shapes.append(dict(
                            type="line",
                            x0=p, x1=p, xref="x",
                            y0=0, y1=1, yref="paper",
                            line=dict(color=color, width=1, dash="solid"),
                            opacity=0.7,
                        ))
        fig.update_layout(shapes=shapes)
        # Legend via invisible scatter traces (one per NlaIII category)
        for cat_lbl, color in seen.items():
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode="lines",
                line=dict(color=color, width=2),
                name=cat_lbl,
            ))

    fig.update_layout(
        title=f"{label} vector coverage \u2014 {contig} ({length:,} bp)",
        xaxis=dict(title=f"{contig} position (bp)", range=[1, length]),
        yaxis=dict(title="per-base depth", range=[0, y_max]),
        width=950, height=620,
        margin=dict(l=60, r=30, t=50, b=50),
        legend=dict(font=dict(size=10)),
    )
    fig.write_image(out_png)
    return (1, length)


def run(sorted_bam, mixed_fa, plasmid_contig, out_prefix):
    """Generate the per-base depth table + linear coverage plot for one contig.

    Returns a dict of the output paths produced.
    """
    out_dir = os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    label = os.path.basename(out_prefix)

    fai = mixed_fa + ".fai"
    length = get_contig_length(fai, plasmid_contig)

    depth_tsv = out_prefix + ".vector_depth.tsv"
    plasmid_depth(sorted_bam, plasmid_contig, depth_tsv)

    png = out_prefix + ".vector_coverage.png"
    if is_nonempty(png):
        log.info(f"Vector coverage plot already present, skipping: {png}")
    else:
        depth = read_depth(depth_tsv, length)
        plot_vector_coverage(depth, length, png, label, plasmid_contig)
        log.info(f"Vector coverage plot -> {png}")

    return {"vector_depth": depth_tsv, "vector_coverage_png": png}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="TLA per-base vector (plasmid) coverage plot from a sorted "
                    "BAM. X-axis spans the full plasmid reference; y-axis linear.")
    parser.add_argument("-b", "--bam", required=True,
                        help="Sorted BAM aligned to the mixed (host+vector) reference")
    parser.add_argument("-r", "--mixed_fa", required=True,
                        help="Mixed reference FASTA (its .fai gives the contig length)")
    parser.add_argument("-c", "--plasmid_contig", required=True,
                        help="Plasmid/vector contig name")
    parser.add_argument("-o", "--out_prefix", required=True,
                        help="Output path prefix (dir is created if needed)")
    args = parser.parse_args()
    paths = run(args.bam, args.mixed_fa, args.plasmid_contig, args.out_prefix)
    for name, path in paths.items():
        log.info(f"  {name}: {path}")
