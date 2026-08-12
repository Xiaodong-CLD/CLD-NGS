#!/usr/bin/env python3

"""TLA-seq Stage 2 -- host-genome coverage (standalone).

From a sorted BAM (mixed host+vector reference) + ``plasmid_contigs.txt``:

    -> host-only coverage bedGraph   (<prefix>.host_only.coverage.bedgraph)
    -> per-chromosome lane plot       (<prefix>.host_coverage.png)
    [-> unplaced-scaffold diagnostic   (<prefix>.signal_contigs.png)]

bedGraph is mapped, non-duplicate, non-secondary alignments (``bedtools
genomecov -bg``); supplementary alignments are kept so split-read segments count
toward the contig they map to, and plasmid contigs are dropped so the trace is
host-only. Binning/plotting come from ``tools/plot_genome_coverage.py``; this
module adds the BAM -> bedGraph step and wires the pieces. Runs standalone (CLI
below) and as Stage 2 of ``pipeline.py``.
"""

import argparse
import os
import subprocess
import sys

from pipeline_io import is_nonempty, bam_ok, commit, get_logger

# The plotting/binning helpers live under tools/ and are self-contained (no
# pipeline_io dependency). Add that dir to the path and import by bare name.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools"))
import plot_genome_coverage as pgc  # noqa: E402

log = get_logger()

# Filter flag for the host coverage track. Supplementary alignments (0x800) are
# KEPT: split-read segments map to different loci, so dropping them would
# undercount coverage at integration junctions. Secondary alignments (0x100)
# are excluded (redundant multi-mapping copies). Duplicate-marked alignments
# (0x400) are KEPT: duplicate removal is not applied here (0x504 would drop
# them). Requires duplicates to be marked upstream (align.py).
_HOST_COVERAGE_FILTER_FLAG = "0x104"


def read_plasmid_contigs(plasmid_contigs_path):
    """Return the set of vector/plasmid contig names (one per line)."""
    with open(plasmid_contigs_path) as fh:
        return {line.strip() for line in fh if line.strip()}


def filter_host_bedgraph_lines(lines, plasmid_set):
    """Yield only the bedGraph lines whose contig is *not* a plasmid contig."""
    for line in lines:
        contig = line.split("\t", 1)[0]
        if contig not in plasmid_set:
            yield line


def host_only_bedgraph(bam, plasmid_contigs_path, out_bedgraph,
                       filter_flag=_HOST_COVERAGE_FILTER_FLAG):
    """Write a host-only coverage bedGraph from ``bam``.

    Pipes ``samtools view -F <filter_flag>`` into ``bedtools genomecov -bg``
    and drops plasmid contigs. No-clobber + atomic: writes a ``.tmp``
    committed only on success.
    """
    if is_nonempty(out_bedgraph):
        log.info(f"Host bedGraph already present, skipping: {out_bedgraph}")
        return out_bedgraph
    if not bam_ok(bam):
        raise FileNotFoundError(f"Sorted BAM missing or truncated: {bam}")

    plasmid_set = read_plasmid_contigs(plasmid_contigs_path)
    tmp = out_bedgraph + ".tmp"
    if os.path.exists(tmp):
        os.remove(tmp)

    log.info(f"Building host-only coverage bedGraph from {os.path.basename(bam)} "
             f"(filter flag {filter_flag}) ...")
    view = subprocess.Popen(
        ["samtools", "view", "-u", "-F", filter_flag, bam],
        stdout=subprocess.PIPE)
    cov = subprocess.Popen(
        ["bedtools", "genomecov", "-bg", "-ibam", "-"],
        stdin=view.stdout, stdout=subprocess.PIPE, text=True)
    view.stdout.close()  # let `view` get SIGPIPE if `cov` exits early

    try:
        with open(tmp, "w") as out:
            for line in filter_host_bedgraph_lines(cov.stdout, plasmid_set):
                out.write(line)
    except BaseException:
        cov.kill()
        view.kill()
        if os.path.exists(tmp):
            os.remove(tmp)
        raise

    cov_rc = cov.wait()
    view_rc = view.wait()
    if view_rc != 0 or cov_rc != 0:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise subprocess.CalledProcessError(
            view_rc or cov_rc, "samtools view | bedtools genomecov")

    if not commit(tmp, out_bedgraph, validator=is_nonempty):
        raise RuntimeError(
            f"Host bedGraph came out empty (no host coverage?): {out_bedgraph}")
    log.info(f"Host-only bedGraph -> {out_bedgraph}")
    return out_bedgraph


def run(bam, fai, plasmid_contigs, out_prefix):
    """Generate the Stage 2 host-coverage outputs for one sample.

    Returns a dict of the output paths produced.
    """
    out_dir = os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    label = os.path.basename(out_prefix)

    bedgraph = out_prefix + ".host_coverage.with_duplicates.bedgraph"
    host_only_bedgraph(bam, plasmid_contigs, bedgraph)

    all_sizes = pgc.read_chrom_sizes(fai, chrom_only=False)
    chrom_map = pgc.read_chrom_map(fai)
    bins = pgc.bin_bedgraph(bedgraph, all_sizes, pgc.BIN_SIZE)

    genome_max = max((arr.max() if len(arr) else 0.0) for arr in bins.values())

    # Main plot: assembled chromosomes only (NC_), in .fai order.
    nc_sizes = {c: L for c, L in all_sizes.items() if c.startswith("NC_")}
    nc_bins = {c: bins[c] for c in nc_sizes}

    line_png = out_prefix + ".host_coverage.png"
    pgc.plot_coverage_tracks(nc_bins, nc_sizes, pgc.BIN_SIZE, line_png, label,
                             chrom_map, global_max=genome_max)

    outputs = {"bedgraph": bedgraph, "host_coverage_png": line_png}

    # Diagnostic: NW_ scaffolds carrying signal.
    signal_sizes = pgc.select_signal_contigs(bins, all_sizes)
    if signal_sizes:
        signal_bins = {c: bins[c] for c in signal_sizes}
        signal_labels = {**chrom_map, **{c: c for c in signal_sizes}}
        signal_png = out_prefix + ".signal_contigs.png"
        pgc.plot_coverage_tracks(signal_bins, signal_sizes, pgc.BIN_SIZE,
                                 signal_png, label + " unplaced scaffolds",
                                 signal_labels, global_max=genome_max)
        outputs["signal_contigs_png"] = signal_png
    else:
        log.info("No unplaced (NW_) scaffolds passed the signal threshold; "
                 "no signal-contigs diagnostic plot written.")

    log.info(f"Stage 2 host coverage complete: {out_prefix}.*")
    return outputs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="TLA Stage 2: host-genome coverage from a sorted BAM "
                    "(host-only bedGraph + whole-genome coverage plot + "
                    "top-bins table).")
    parser.add_argument("-b", "--bam", required=True,
                        help="Sorted BAM aligned to the mixed (host+vector) reference")
    parser.add_argument("-p", "--plasmid_contigs", required=True,
                        help="plasmid_contigs.txt (vector contig names; from build_reference)")
    parser.add_argument("-f", "--fai", required=True,
                        help="Host genome .fai (chromosome sizes/order; chrom_map.tsv "
                             "is auto-loaded if beside it)")
    parser.add_argument("-o", "--out_prefix", required=True,
                        help="Output path prefix (dir is created if needed)")
    args = parser.parse_args()
    paths = run(args.bam, args.fai, args.plasmid_contigs, args.out_prefix)
    for name, path in paths.items():
        log.info(f"  {name}: {path}")
