#!/usr/bin/env python3

"""Split reportable GRIDSS host-vector junctions into resolved and unresolved TSVs.

The upstream ``integration_sites.tsv`` remains the complete host-vector evidence
file. Junctions classified as ``both_ends_near_nlaiii`` are omitted from the two
user-facing reports but remain preserved upstream.

Outputs
-------
<prefix>.resolved_integration_sites.tsv
    One row per reportable host-vector junction in a site cluster containing
    both host-facing arms.

<prefix>.unresolved_host_vector_junctions.tsv
    One row per remaining reportable host-vector junction. One-end-near-NlaIII
    and neither-end-near-NlaIII junctions are retained with their annotation.

<prefix>.integration_report_summary.tsv
    Count summary for upstream junction classes, NlaIII categories, and final
    resolved/unresolved site and junction totals.

``site_id`` is retained only to link junctions assigned to the same
local candidate site cluster. Resolution is based on whether both host-facing
arms remain after NlaIII filtering.
"""

import argparse
import csv
import os
import sys

from pipeline_io import get_logger

log = get_logger()

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
HOST_WINDOW = 1000        # bp; single-linkage site clustering distance
_NLAIII_COLUMN = "nlaiii_category"
_NLAIII_BOTH = "both_ends_near_nlaiii"

# Columns present in the combined gridss_breakpoints.tsv that are not part of
# the host-vector integration schema and must be excluded when reading from it.
_COMBINED_ONLY_COLS = frozenset({
    "breakpoint_type", "contig_a", "pos_a", "contig_b", "pos_b",
    "host_a", "host_b", "contig", "pos", "is_vector",
    "nlaiii_distance_1", "nlaiii_distance_2",
})


# ---------------------------------------------------------------------------
# Shared helpers (inlined from gridss_report)
# ---------------------------------------------------------------------------

def _to_int(value, default=0):
    try:
        return int(value)
    except (TypeError, ValueError):
        try:
            return int(float(value))
        except (TypeError, ValueError):
            return default


def _host_arm(breakend):
    """'right' if the VCF ALT starts with a bracket, else 'left'."""
    return "right" if breakend and breakend[0] in "[]" else "left"


def read_integration_sites(path):
    """Read host-vector rows from a gridss_breakpoints.tsv or legacy integration_sites.tsv.

    When the file contains a ``breakpoint_type`` column (combined format), only
    rows with ``breakpoint_type == "integration_sites"`` are returned, projected
    to the integration-specific columns. Legacy single-type files are read as-is.
    Also returns whether the rows have a ``nlaiii_category`` column.
    """
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        all_rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if "breakpoint_type" in fieldnames:
        all_rows = [r for r in all_rows if r.get("breakpoint_type") == "integration_sites"]
        fieldnames = [col for col in fieldnames if col not in _COMBINED_ONLY_COLS]
        all_rows = [{col: row.get(col, "") for col in fieldnames} for row in all_rows]
    has_nlaiii = _NLAIII_COLUMN in fieldnames
    return all_rows, has_nlaiii


def drop_both_ends_near_nlaiii(rows):
    """Drop both-end NlaIII junctions (likely religation artifacts). Returns (kept, n_dropped)."""
    kept = [r for r in rows if r.get(_NLAIII_COLUMN) != _NLAIII_BOTH]
    return kept, len(rows) - len(kept)


def _cluster(rows):
    """Single-linkage cluster junctions into sites by host position (HOST_WINDOW bp)."""
    ordered = sorted(rows, key=lambda r: (r["host_contig"], _to_int(r["host_pos"])))
    sites = []
    prev_contig = prev_pos = None
    for row in ordered:
        contig, pos = row["host_contig"], _to_int(row["host_pos"])
        if not sites or prev_contig != contig or (pos - prev_pos) > HOST_WINDOW:
            sites.append({"contig": contig, "rows": []})
        sites[-1]["rows"].append(row)
        prev_contig, prev_pos = contig, pos
    for idx, site in enumerate(sites, start=1):
        site["site_id"] = f"site_{idx:03d}"
    return sites

_META_COLUMNS = ("site_id", "arm")


def _build_sites(rows, has_nlaiii):
    """Cluster reportable junctions and annotate host-facing arms.

    Arm orientation is read directly from each raw row's ``breakend`` field
    using ``_host_arm``, preserving the existing orientation convention.
    """
    sites = _cluster(rows)

    for site in sites:
        has_left = any(_host_arm(r.get("breakend", "")) == "left" for r in site["rows"])
        has_right = any(_host_arm(r.get("breakend", "")) == "right" for r in site["rows"])
        site["_resolved"] = has_left and has_right
        site["_min_pos"] = min(int(row["host_pos"]) for row in site["rows"])

    sites.sort(key=lambda site: (site["contig"], site["_min_pos"], site["site_id"]))
    return sites


def _records_for_status(sites, row_columns, want_resolved):
    """Yield every junction row from sites matching ``want_resolved`` exactly once."""
    for site in sites:
        if site["_resolved"] != want_resolved:
            continue
        for row in site["rows"]:
            record = {
                "site_id": site["site_id"],
                "arm": _host_arm(row.get("breakend", "")),
            }
            record.update((column, row.get(column, "")) for column in row_columns)
            yield record


def _write_tsv(path, records, row_columns):
    """Atomically write a single-header junction TSV and return its row count."""
    columns = list(_META_COLUMNS) + list(row_columns)
    tmp = path + ".tmp"
    n_rows = 0
    with open(tmp, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for record in records:
            fh.write("\t".join(str(record.get(column, "")) for column in columns) + "\n")
            n_rows += 1
    os.replace(tmp, path)
    return n_rows


def run(integration_sites_tsv, out_prefix):
    """Write resolved, unresolved, and count-summary integration reports."""
    resolved_path = f"{out_prefix}.resolved_sites.tsv"
    unresolved_path = f"{out_prefix}.unresolved_junctions.tsv"

    all_rows, has_nlaiii = read_integration_sites(integration_sites_tsv)
    reportable_rows = all_rows
    n_dropped = 0
    if has_nlaiii:
        reportable_rows, n_dropped = drop_both_ends_near_nlaiii(all_rows)
        if n_dropped:
            log.info(
                "Excluded %d both-end-NlaIII junctions from user-facing reports; "
                "records remain in %s",
                n_dropped,
                integration_sites_tsv,
            )

    sites = _build_sites(reportable_rows, has_nlaiii)
    # Use the actual input fieldnames so the header matches the data rows.
    row_columns = tuple(all_rows[0].keys()) if all_rows else ()

    n_resolved = _write_tsv(
        resolved_path,
        _records_for_status(sites, row_columns, want_resolved=True),
        row_columns,
    )
    n_unresolved = _write_tsv(
        unresolved_path,
        _records_for_status(sites, row_columns, want_resolved=False),
        row_columns,
    )

    if n_resolved + n_unresolved != len(reportable_rows):
        raise RuntimeError(
            "Resolved/unresolved reports did not partition all reportable junctions: "
            f"{n_resolved} + {n_unresolved} != {len(reportable_rows)}"
        )

    log.info(
        "Integration reports: %d resolved junction rows, %d unresolved junction "
        "rows, %d excluded both-end-NlaIII rows",
        n_resolved,
        n_unresolved,
        n_dropped,
    )
    log.info("  resolved -> %s", resolved_path)
    log.info("  unresolved -> %s", unresolved_path)

    return {
        "resolved": resolved_path,
        "unresolved": unresolved_path,
        "counts": {
            "n_resolved_sites": n_resolved,
            "n_unresolved_junctions": n_unresolved,
            "excluded_both_ends_near_nlaiii": n_dropped,
        },
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Split GRIDSS host-vector junctions into resolved and unresolved TSVs."
    )
    parser.add_argument(
        "-i",
        "--integration-sites",
        required=True,
        help="Upstream GRIDSS host-vector junction table (*.integration_sites.tsv).",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        required=True,
        help="Output prefix for resolved, unresolved, and summary TSVs.",
    )
    args = parser.parse_args(argv)
    run(args.integration_sites, args.out_prefix)


if __name__ == "__main__":
    main()