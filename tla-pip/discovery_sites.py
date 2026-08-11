#!/usr/bin/env python3
"""Coarse host-enrichment discovery regions from a host bedGraph + GRIDSS outputs (TLA).

Bins the bedGraph, merges enriched windows into candidate regions, matches
each region against resolved and unresolved GRIDSS outputs, and writes one
<prefix>.discovery_sites.tsv.

Runs as a step in pipeline.py after host_coverage and gridss_split_report.
"""

import csv
import os
from collections import defaultdict

from pipeline_io import get_logger

log = get_logger()

DISCOVERY_COVERAGE_CUTOFF = 5
DISCOVERY_BIN_SIZE_BP = 10000000
DISCOVERY_MAX_GAP_BINS = 1

_HEADER = (
    "candidate_region_id", "contig", "start", "end", "candidate_region",
    "max_depth", "covered_bp_above_cutoff",
    "n_gridss_breakends", "n_resolved_sites", "n_unresolved_junctions",
    "resolved_site_ids", "unresolved_junction_ids", "support_class",
)


def _flush(regions, contig, b_start, b_end, bin_size, bin_max, bin_covered):
    r_start = b_start * bin_size
    r_end = (b_end + 1) * bin_size
    r_max = max(bin_max[(contig, b)] for b in range(b_start, b_end + 1))
    r_cov = sum(bin_covered[(contig, b)] for b in range(b_start, b_end + 1))
    regions.append({
        "contig": contig, "start": r_start, "end": r_end,
        "max_depth": r_max, "covered_bp_above_cutoff": r_cov,
    })


def _build_regions(bedgraph_path, bin_size, cutoff, max_gap):
    """Single bedGraph pass: find active bins, merge, compute stats.

    Returns list of dicts: contig, start, end (0-based half-open),
    max_depth, covered_bp_above_cutoff. Active bins are those containing at
    least one bedGraph interval with depth > cutoff. Adjacent active bins and
    bins separated by at most max_gap inactive bins are merged into one region.
    """
    bin_max = defaultdict(float)
    bin_covered = defaultdict(int)

    with open(bedgraph_path) as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            contig = parts[0]
            s, e, v = int(parts[1]), int(parts[2]), float(parts[3])
            if v <= cutoff:
                continue
            b0, b1 = s // bin_size, (e - 1) // bin_size
            for b in range(b0, b1 + 1):
                lo = max(s, b * bin_size)
                hi = min(e, (b + 1) * bin_size)
                key = (contig, b)
                if v > bin_max[key]:
                    bin_max[key] = v
                bin_covered[key] += hi - lo

    contig_bins = defaultdict(list)
    for contig, b in bin_max:
        contig_bins[contig].append(b)

    regions = []
    for contig in sorted(contig_bins):
        sorted_b = sorted(contig_bins[contig])
        block_s = sorted_b[0]
        block_e = sorted_b[0]
        for b in sorted_b[1:]:
            if b - block_e <= max_gap + 1:
                block_e = b
            else:
                _flush(regions, contig, block_s, block_e, bin_size, bin_max, bin_covered)
                block_s = block_e = b
        _flush(regions, contig, block_s, block_e, bin_size, bin_max, bin_covered)

    return regions


def _read_tsv(path):
    if not os.path.exists(path):
        log.warning("GRIDSS output not found, skipping match: %s", path)
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _annotate(regions, resolved_rows, unresolved_rows):
    """Attach GRIDSS evidence to each region.

    Coordinate note: host_pos in GRIDSS outputs is 1-based (VCF convention);
    bedGraph regions are 0-based half-open [start, end).
    A 1-based pos p is inside [start, end) iff start < p <= end.
    """
    n = len(regions)
    res_ids = [set() for _ in range(n)]   # unique site_ids from resolved rows
    unres_ids = [[] for _ in range(n)]     # one entry per unresolved row

    for row in resolved_rows:
        contig = row.get("host_contig", "")
        try:
            pos = int(row.get("host_pos", 0))
        except (ValueError, TypeError):
            continue
        for i, reg in enumerate(regions):
            if reg["contig"] == contig and reg["start"] < pos <= reg["end"]:
                res_ids[i].add(row.get("site_id", ""))

    for row in unresolved_rows:
        contig = row.get("host_contig", "")
        try:
            pos = int(row.get("host_pos", 0))
        except (ValueError, TypeError):
            continue
        for i, reg in enumerate(regions):
            if reg["contig"] == contig and reg["start"] < pos <= reg["end"]:
                unres_ids[i].append(row.get("site_id", ""))

    for i, reg in enumerate(regions):
        rs = sorted(res_ids[i])
        ur = unres_ids[i]
        n_res = len(rs)
        n_unres = len(ur)
        reg["n_resolved_sites"] = n_res
        reg["n_unresolved_junctions"] = n_unres
        reg["n_gridss_breakends"] = n_res + n_unres
        reg["resolved_site_ids"] = ",".join(rs) if rs else "."
        reg["unresolved_junction_ids"] = ",".join(ur) if ur else "."
        if n_res > 0:
            reg["support_class"] = "resolved_supported"
        elif n_unres > 0:
            reg["support_class"] = "unresolved_supported"
        else:
            reg["support_class"] = "coverage_only"

    return regions


def run(bedgraph, resolved_tsv, unresolved_tsv, out_prefix):
    """Generate <out_prefix>.discovery_sites.tsv.

    Returns count dict: n_discovery_regions, n_gridss_supported_discovery_regions,
    n_resolved_discovery_regions.
    """
    out_path = out_prefix + ".discovery_sites.tsv"

    regions = _build_regions(bedgraph, DISCOVERY_BIN_SIZE_BP,
                             DISCOVERY_COVERAGE_CUTOFF, DISCOVERY_MAX_GAP_BINS)
    log.info(
        "Coverage enrichment: %d coarse regions (cutoff=%s, bin=%dbp, max_gap=%d bins)",
        len(regions), DISCOVERY_COVERAGE_CUTOFF, DISCOVERY_BIN_SIZE_BP, DISCOVERY_MAX_GAP_BINS,
    )

    resolved_rows = _read_tsv(resolved_tsv)
    unresolved_rows = _read_tsv(unresolved_tsv)
    regions = _annotate(regions, resolved_rows, unresolved_rows)

    tmp = out_path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(_HEADER) + "\n")
        for idx, reg in enumerate(regions, start=1):
            region_id = f"DISCOVERY_{idx:04d}"
            fh.write("\t".join(str(x) for x in [
                region_id,
                reg["contig"],
                reg["start"],
                reg["end"],
                f"{reg['contig']}:{reg['start']}-{reg['end']}",
                f"{reg['max_depth']:.2f}",
                reg["covered_bp_above_cutoff"],
                reg["n_gridss_breakends"],
                reg["n_resolved_sites"],
                reg["n_unresolved_junctions"],
                reg["resolved_site_ids"],
                reg["unresolved_junction_ids"],
                reg["support_class"],
            ]) + "\n")
    os.replace(tmp, out_path)
    log.info("Discovery sites -> %s", out_path)

    n_total = len(regions)
    n_gridss = sum(1 for r in regions if r["n_gridss_breakends"] > 0)
    n_resolved = sum(1 for r in regions if r["n_resolved_sites"] > 0)
    log.info("  %d total, %d GRIDSS-supported, %d resolved-supported",
             n_total, n_gridss, n_resolved)

    return {
        "n_discovery_regions": n_total,
        "n_gridss_supported_discovery_regions": n_gridss,
        "n_resolved_discovery_regions": n_resolved,
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="TLA discovery sites: coarse coverage regions + GRIDSS overlap."
    )
    parser.add_argument("-b", "--bedgraph", required=True,
                        help="Host bedGraph (output of host_coverage.py)")
    parser.add_argument("-r", "--resolved", required=True,
                        help="<prefix>.resolved_integration_sites.tsv")
    parser.add_argument("-u", "--unresolved", required=True,
                        help="<prefix>.unresolved_host_vector_junctions.tsv")
    parser.add_argument("-o", "--out_prefix", required=True,
                        help="Output path prefix")
    args = parser.parse_args()
    run(args.bedgraph, args.resolved, args.unresolved, args.out_prefix)
