#!/usr/bin/env python3
"""Raw-depth coverage plot for one simulated condition (diagnostic).

Reads the per-base raw depth produced by the pipeline's samtools step
(`samtools depth -a` -> <sample>.depth.txt) and renders depth vs position as a
PNG, marking the injected indel position and the 5'/3' read-covered zones. With
a small primer-to-primer amplicon and 2x75 reads the profile is expected to
peak at both ends and dip in the center insert gap; this plot makes that
visible for diagnosis.

Standalone (does NOT reuse plot_coverage_plots.py, which needs a VarScan
consensus that this benchmark skips).
"""

import argparse
import os
import sys

import plotly.graph_objects as go

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config  # noqa: E402


def read_depth(depth_path):
    """Return (positions, depths) from a samtools depth -a file."""
    positions, depths = [], []
    with open(depth_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            positions.append(int(parts[1]))
            depths.append(int(parts[2]))
    return positions, depths


def plot_depth(depth_path, out_png, title=None, indel_pos=None,
               read_len=config.READ_LEN, amplicon_len=None):
    positions, depths = read_depth(depth_path)
    if not positions:
        print(f"[coverage_plot] No depth data in {depth_path}, skipping")
        return None

    # Fall back to the observed max position when amplicon_len is not supplied.
    if amplicon_len is None:
        amplicon_len = max(positions)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=positions, y=depths, mode="lines",
                             name="raw depth", line=dict(color="#1f77b4")))

    # Read-covered zones (5' and 3' ends of the amplicon).
    fig.add_vrect(x0=1, x1=read_len, fillcolor="green", opacity=0.07,
                  line_width=0, annotation_text="5' reads",
                  annotation_position="top left")
    fig.add_vrect(x0=amplicon_len - read_len + 1, x1=amplicon_len,
                  fillcolor="green", opacity=0.07, line_width=0,
                  annotation_text="3' reads", annotation_position="top right")

    # Injected indel position (skipped for WT-only controls).
    if indel_pos is not None:
        fig.add_vline(x=indel_pos, line=dict(color="red", dash="dash"),
                      annotation_text="indel", annotation_position="top")

    fig.update_layout(
        title=title or os.path.basename(depth_path),
        xaxis_title="amplicon position (bp)",
        yaxis_title="raw depth",
        template="plotly_white",
        width=900, height=400,
        showlegend=False,
    )

    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    fig.write_image(out_png)
    return out_png


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("depth_file", help="Path to <sample>.depth.txt")
    parser.add_argument("-o", "--out", required=True, help="Output PNG path")
    parser.add_argument("-t", "--title", default=None)
    parser.add_argument("--indel-pos", type=int, default=config.INDEL_POS)
    args = parser.parse_args()

    out = plot_depth(args.depth_file, args.out, title=args.title,
                     indel_pos=args.indel_pos)
    if out:
        print(f"[coverage_plot] Wrote {out}")


if __name__ == "__main__":
    main()
