#!/usr/bin/env python3

"""Binned host-coverage plot + top-bins table from a host-only bedGraph (TLA QC).

A raw per-base bedGraph is too large for IGV, so it is binned into fixed windows
and drawn as a stacked per-chromosome lane plot; a TLA integration locus shows
as a tall coverage spike. 

Inputs: a 4-column bedGraph (``chrom start end depth``)
and the host ``.fai``. 

Outputs: ``<out_prefix>.host_coverage.png`` and
``<out_prefix>.top_bins.tsv``. A ``chrom_map.tsv`` beside the ``.fai`` (if
present) maps accessions to chromosome labels (e.g. ``chr4``).

Fixed settings are module constants.
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless: write PNG without a display
import matplotlib.pyplot as plt

BIN_SIZE = 100_000      # bin width in bp
TOP_N = 50              # number of top bins written to the TSV
CHROM_ONLY = True       # only assembled chromosomes (NC_ accessions)
CHROM_MAP_NAME = "chrom_map.tsv"  # looked for next to the .fai

# Signal-contigs diagnostic: chr1 in CriGri-PICRH-1.0 exists only as unlocalized
# NW_ scaffolds (excluded from the NC_-only main plot), so an integration there
# would be invisible. A NW_ scaffold is shown if its peak bin >= SIGNAL_MIN_DEPTH
# (mean depth per 100kb bin); capped at SIGNAL_MAX_TRACKS lanes (strongest
# first). An absolute floor is used rather than a fraction of the genome-wide
# peak: the NW_ peak distribution is heavily background (lib1: 90% of scaffolds
# peak < 0.7), and a fractional cutoff collapses to noise whenever the genome
# peak is small. A depth floor cleanly separates real pile-up from background
# regardless of how tall the main locus is. The genome-wide thin background
# (~1-2x per bin) is expected for TLA proximity ligation (3D contacts, genomic
# carryover, repeat/multi-mapping), not spurious coverage to remove.
# TODO(investigate later): a per-library relative cutoff (e.g. a multiple of the
# NW_ background p99, or a Poisson background model) would adapt the floor to
# library depth automatically instead of the fixed 5.0.
UNPLACED_PREFIX = "NW_"
SIGNAL_MIN_DEPTH = 5.0
SIGNAL_MAX_TRACKS = 50


def read_chrom_sizes(fai_path, chrom_only):
    """Ordered ``{name: length}`` from a .fai. In chromosome-only mode, keep
    only assembled chromosomes (``NC_``) and drop ``NW_`` scaffolds."""
    sizes = {}
    with open(fai_path) as fh:
        for line in fh:
            fields = line.split("\t")
            name, length = fields[0], int(fields[1])
            if chrom_only and not name.startswith("NC_"):
                continue
            sizes[name] = length
    if not sizes:
        raise ValueError(
            f"No sequences selected from {fai_path} "
            f"(chrom_only={chrom_only}). Check the .fai naming.")
    return sizes


def read_chrom_map(fai_path):
    """Return ``{accession: label}`` from a ``chrom_map.tsv`` (header +
    accession/chromosome rows) beside the ``.fai``, or ``{}`` if absent."""
    map_path = os.path.join(os.path.dirname(fai_path), CHROM_MAP_NAME)
    mapping = {}
    if not os.path.exists(map_path):
        return mapping
    with open(map_path) as fh:
        next(fh, None)  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0]:
                mapping[parts[0]] = parts[1]
    return mapping


def chrom_label(accession, chrom_map):
    """Mapped chromosome name, else the trailing accession digits."""
    if accession in chrom_map:
        return chrom_map[accession]
    return accession.split(".")[0].split("_")[-1]


def bin_bedgraph(bedgraph_path, sizes, bin_size):
    """Bin per-base depth into fixed windows -> ``{chrom: mean-depth array}``.

    Each interval adds ``value * overlap_bp`` to every bin it touches; the
    per-bin sum is divided by the bin width to give mean depth. Intervals on
    contigs without a known size are skipped.
    """
    bins = {c: np.zeros(L // bin_size + 1) for c, L in sizes.items()}
    with open(bedgraph_path) as fh:
        for line in fh:
            c, s, e, v = line.split("\t")
            arr = bins.get(c)
            if arr is None:
                continue
            s, e, v = int(s), int(e), float(v)
            b0, b1 = s // bin_size, (e - 1) // bin_size
            for b in range(b0, b1 + 1):
                lo = max(s, b * bin_size)
                hi = min(e, (b + 1) * bin_size)
                arr[b] += v * (hi - lo)
    for c in bins:
        bins[c] /= bin_size
    return bins


def write_top_bins(bins, bin_size, out_tsv, top_n, chrom_map):
    """Write the top-N bins by mean depth to a TSV
    (chromosome/accession/start/end/mean_depth)."""
    rows = []
    for c, arr in bins.items():
        for b, depth in enumerate(arr):
            rows.append((c, b * bin_size, (b + 1) * bin_size, float(depth)))
    rows.sort(key=lambda r: r[3], reverse=True)
    with open(out_tsv, "w") as fh:
        fh.write("chromosome\taccession\tstart\tend\tmean_depth\n")
        for c, start, end, depth in rows[:top_n]:
            fh.write(f"{chrom_label(c, chrom_map)}\t{c}\t{start}\t{end}\t{depth:.2f}\n")


def plot_coverage_tracks(bins, sizes, bin_size, out_png, sample_label, chrom_map,
                         track_fraction=0.85, global_max=None, bins_with_dup=None):
    """Draw a stacked per-chromosome lane coverage plot (Cegentis Fig-8 style).

    Each contig is a horizontal baseline drawn to its true length on a shared
    x-axis (Mb), with coverage as vertical spikes; an integration locus shows
    as a tall spike. Contig order sets top-to-bottom (first = top lane). All
    lanes share one vertical scale (the figure's tallest bin) so peak heights
    are comparable; the tallest peak fills a configurable fraction of a lane.

    Pass ``global_max`` (the genome-wide peak depth) to anchor that scale to a
    value shared across figures, so the main plot and the signal-contigs
    diagnostic use identical depth-per-pixel. Defaults to this figure's own
    tallest bin.

    Pass ``bins_with_dup`` to overlay a second track (duplicate-retained
    alignments) in steelblue alongside ``bins`` (duplicate-excluded) in
    crimson; both tracks use the same scale.
    """
    order = list(sizes.keys())
    n = len(order)

    # Shared vertical scale: a peak height means the same depth on every lane.
    if global_max is None:
        global_max = max((arr.max() if len(arr) else 0.0) for arr in bins.values())
    scale = track_fraction / (global_max or 1.0)

    overlay = bins_with_dup is not None
    color_nd = "crimson" if overlay else "black"
    lw = 0.5 if overlay else 0.4
    alpha = 0.7 if overlay else 1.0

    fig_h = max(3.0, 0.32 * n + 1.0)
    plt.figure(figsize=(14, fig_h))
    yticks, ylabels = [], []
    # Reverse so the first contig in ``sizes`` lands on the top lane (row 0 is
    # drawn at the bottom).
    for r, c in enumerate(reversed(order)):
        length_mb = sizes[c] / 1e6
        plt.hlines(r, 0, length_mb, color="0.6", lw=0.6)  # chromosome baseline
        arr = bins[c]
        if len(arr):
            xs = (np.arange(len(arr)) + 0.5) * bin_size / 1e6
            heights = arr * scale
            nz = heights > 0  # only draw covered bins (keeps the figure light)
            if nz.any():
                plt.vlines(xs[nz], r, r + heights[nz], color=color_nd, lw=lw, alpha=alpha)
        if overlay:
            arr_wd = bins_with_dup.get(c)
            if arr_wd is not None and len(arr_wd):
                xs_wd = (np.arange(len(arr_wd)) + 0.5) * bin_size / 1e6
                heights_wd = arr_wd * scale
                nz_wd = heights_wd > 0
                if nz_wd.any():
                    plt.vlines(xs_wd[nz_wd], r, r + heights_wd[nz_wd],
                               color="steelblue", lw=lw, alpha=alpha)
        yticks.append(r)
        ylabels.append(chrom_label(c, chrom_map))

    plt.yticks(yticks, ylabels, fontsize=7)
    plt.ylim(-0.5, n - 0.5 + track_fraction)
    plt.xlim(left=0)
    plt.xlabel("chromosomal position (Mb)")
    if overlay:
        from matplotlib.lines import Line2D
        plt.legend(
            handles=[
                Line2D([0], [0], color="steelblue", alpha=alpha, lw=1.5,
                       label="with duplicates"),
                Line2D([0], [0], color=color_nd, alpha=alpha, lw=1.5,
                       label="without duplicates"),
            ],
            loc="upper right", fontsize=7,
        )
        plt.title(f"{sample_label} -- host coverage with vs. without "
                  f"duplicate-marked reads")
    else:
        plt.title(f"{sample_label} host coverage -- peaks indicate candidate "
                  f"integration loci")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


def select_signal_contigs(bins, sizes, min_depth=SIGNAL_MIN_DEPTH,
                          max_tracks=SIGNAL_MAX_TRACKS, prefix=UNPLACED_PREFIX):
    """Pick unplaced (NW_) scaffolds with a coverage peak worth plotting.

    A scaffold qualifies when its tallest bin reaches at least ``min_depth``
    (mean depth per bin, default SIGNAL_MIN_DEPTH). An absolute floor is used
    instead of a fraction of the genome-wide peak because the NW_ peak
    distribution is dominated by background, and a fractional cutoff collapses
    to noise when the genome peak is small. Qualifying scaffolds are returned
    strongest-first and capped at the configured track limit (default
    SIGNAL_MAX_TRACKS). Returns ``{contig: length_bp}`` ordered for plotting, or
    ``{}`` if nothing clears the bar.
    """
    hits = []
    for c in sizes:
        if not c.startswith(prefix):
            continue
        arr = bins.get(c)
        peak = arr.max() if (arr is not None and len(arr)) else 0.0
        if peak >= min_depth:
            hits.append((c, peak))
    hits.sort(key=lambda t: t[1], reverse=True)
    return {c: sizes[c] for c, _ in hits[:max_tracks]}


def main():
    parser = argparse.ArgumentParser(
        description="Genome-wide binned coverage plot + top-bins table from a "
                    "host-only bedGraph (TLA integration-site QC).")
    parser.add_argument("-b", "--bedgraph", required=True,
                        help="Host-only bedGraph (chrom start end depth)")
    parser.add_argument("-f", "--fai", required=True,
                        help="Host genome .fai (chromosome sizes/order)")
    parser.add_argument("-o", "--out_prefix", required=True,
                        help="Output path prefix (dir is created if needed)")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    label = os.path.basename(args.out_prefix)

    sizes = read_chrom_sizes(args.fai, chrom_only=CHROM_ONLY)
    chrom_map = read_chrom_map(args.fai)
    bins = bin_bedgraph(args.bedgraph, sizes, BIN_SIZE)

    out_png = f"{args.out_prefix}.host_coverage.png"
    out_tsv = f"{args.out_prefix}.top_bins.tsv"
    plot_coverage_tracks(bins, sizes, BIN_SIZE, out_png, label, chrom_map)
    write_top_bins(bins, BIN_SIZE, out_tsv, TOP_N, chrom_map)

    print(f"wrote {out_png}")
    print(f"wrote {out_tsv}")


if __name__ == "__main__":
    main()
