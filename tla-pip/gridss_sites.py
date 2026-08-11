#!/usr/bin/env python3

"""Parse GRIDSS2 PASS breakends into four classified tables.

Reads PASS BND records from the GRIDSS VCF and, in a single pass, classifies
every breakend junction by the contigs of its two anchors (the record ``CHROM``
and the mate contig from the ALT bracket), using the vector/plasmid contig list
(everything else is host, both ``NC_`` assembled chromosomes and ``NW_`` unplaced
scaffolds). It writes four TSVs:

    integration_sites.tsv   host<->vector junctions
    vector_sv.tsv           vector<->vector (within-vector SV candidates)
    host_host_bnd.tsv       host<->host (kept for manual inspection)
    single_breakend.tsv     single-ended breakends

Each two-sided junction appears twice in the VCF (a breakend + its ``MATEID``
partner sharing one ``EVENT``); rows are grouped by ``EVENT`` so each junction is
emitted once, host<->vector anchored on the host side. The host breakpoint keeps
its ``CIPOS`` interval (the microhomology window, its width equals ``HOMLEN``)
rather than a single base. The raw ALT string is carried through as the
``breakend`` column (VCF BND notation, e.g. ``T]D08_HI_HS_AvB1:7626]``). The VCF
is BGZF (gzip-readable), so it is parsed as plain text.
"""

import argparse
import gzip
import os
import re

from pipeline_io import get_logger, is_nonempty, write_tsv
import restriction_site_annotation as rsa

log = get_logger()

# GRIDSS BND ALT: a single mate bracket pair, e.g. ``T[chr:123[`` / ``]chr:123]A``.
_ALT_BRACKET_RE = re.compile(r"([\[\]])([^\[\]]+)\1")


def read_plasmid_contigs(plasmid_contigs_path):
    """Return the set of vector/plasmid contig names (one per line)."""
    with open(plasmid_contigs_path) as fh:
        return {line.strip() for line in fh if line.strip()}


def _count_tsv_rows(path):
    """Return the number of data rows in a TSV (total lines minus the header)."""
    with open(path) as fh:
        n_lines = sum(1 for _ in fh)
    return max(n_lines - 1, 0)


def _parse_info(info_field):
    """Parse a VCF INFO string into a dict; bare flags map to ``True``."""
    info = {}
    for entry in info_field.split(";"):
        if not entry:
            continue
        if "=" in entry:
            key, val = entry.split("=", 1)
            info[key] = val
        else:
            info[entry] = True
    return info


def _num(value, cast, default=""):
    """Cast ``value`` with ``cast``; return ``default`` on missing/bad input."""
    if value in (None, "", True):
        return default
    try:
        return cast(value)
    except (TypeError, ValueError):
        return default


def _parse_alt(alt):
    """Parse the mate anchor out of a GRIDSS BND ALT string.

    Returns ``(mate_contig, mate_pos)`` for a paired breakend, or ``(None,
    None)`` for a single breakend (ALT carries a ``.`` and no mate bracket). The
    full orientation (which arm is kept, whether the mate is reverse-
    complemented) is preserved by carrying the raw ALT string through as the
    ``breakend`` column rather than re-encoding it here.
    """
    match = _ALT_BRACKET_RE.search(alt)
    if match is None:
        return None, None
    mate_contig, mate_pos = match.group(2).rsplit(":", 1)
    return mate_contig, int(mate_pos)


def _parse_record(cols):
    """Parse one PASS VCF data line (already split on tabs) into a dict.

    Pulls the coordinate, ALT-derived mate + orientation, and the INFO support
    fields reported. Returns ``None`` if the record is not a usable BND line.
    """
    chrom, pos, _id, _ref, alt, qual = cols[0], cols[1], cols[2], cols[3], cols[4], cols[5]
    info = _parse_info(cols[7])
    mate_contig, mate_pos = _parse_alt(alt)
    beid = info.get("BEID", "")
    return {
        "chrom": chrom,
        "pos": int(pos),
        "breakend": alt,
        "qual": _num(qual, float),
        "mate_contig": mate_contig,
        "mate_pos": mate_pos,
        "cipos": info.get("CIPOS", "0,0"),
        "homlen": _num(info.get("HOMLEN"), int),
        "homseq": info.get("HOMSEQ", ""),
        "vf": _num(info.get("VF"), int),
        "sr": _num(info.get("SR"), int),
        "rp": _num(info.get("RP"), int),
        "assr": _num(info.get("ASSR"), int),
        "asrp": _num(info.get("ASRP"), int),
        "n_beid": len(beid.split(",")) if beid else 0,
        "imprecise": "IMPRECISE" in info,
        "sb": _num(info.get("SB"), float),
        "event": info.get("EVENT", ""),
        "is_single": mate_contig is None,
    }


def _read_pass_records(vcf_path):
    """Yield parsed PASS BND records from a (BGZF) GRIDSS VCF."""
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8 or cols[6] != "PASS":
                continue
            record = _parse_record(cols)
            if record is not None:
                yield record


def _classify(records, vector_set):
    """Group PASS records by EVENT and split them into the four tables.

    Returns ``(integration, vector_sv, host_host, single)`` lists of row tuples,
    each already ordered for its table's header. Single breakends are emitted per
    record (deduped by contig+pos); two-sided junctions are emitted once per
    EVENT, host<->vector rows anchored on the host side.
    """
    by_event = {}
    singles = {}
    for rec in records:
        if rec["is_single"]:
            singles.setdefault((rec["chrom"], rec["pos"]), rec)
        else:
            by_event.setdefault(rec["event"], []).append(rec)

    integration, vector_sv, host_host = [], [], []
    for rows in by_event.values():
        rep = rows[0]
        c1_vector = rep["chrom"] in vector_set
        c2_vector = rep["mate_contig"] in vector_set

        if c1_vector and c2_vector:
            vector_sv.append((
                rep["chrom"], rep["pos"], rep["mate_contig"], rep["mate_pos"],
                rep["breakend"], rep["qual"], rep["vf"], rep["sr"], rep["rp"],
                rep["assr"], rep["asrp"], rep["homlen"], rep["homseq"],
                rep["n_beid"], rep["imprecise"], rep["sb"], rep["event"]))
        elif not c1_vector and not c2_vector:
            host_host.append((
                rep["chrom"], rep["pos"], rep["cipos"],
                rep["mate_contig"], rep["mate_pos"], rep["breakend"],
                rep["qual"], rep["vf"], rep["sr"], rep["rp"], rep["assr"],
                rep["asrp"], rep["homlen"], rep["homseq"], rep["n_beid"],
                rep["imprecise"], rep["sb"], rep["event"]))
        else:
            # host<->vector: anchor on the host side. Prefer the row whose CHROM
            # is host (its CIPOS is the host-breakpoint interval); fall back to
            # the vector-CHROM row and read the host anchor from its mate.
            host_row = next((r for r in rows if r["chrom"] not in vector_set), None)
            if host_row is not None:
                host_contig, host_pos = host_row["chrom"], host_row["pos"]
                vector_contig, vector_pos = host_row["mate_contig"], host_row["mate_pos"]
                anchor = host_row
            else:
                host_contig, host_pos = rep["mate_contig"], rep["mate_pos"]
                vector_contig, vector_pos = rep["chrom"], rep["pos"]
                anchor = rep
            integration.append((
                host_contig, host_pos, anchor["cipos"],
                vector_contig, vector_pos, anchor["breakend"], anchor["qual"],
                anchor["vf"], anchor["sr"], anchor["rp"], anchor["assr"],
                anchor["asrp"], anchor["homlen"], anchor["homseq"],
                anchor["n_beid"], anchor["imprecise"], anchor["sb"],
                anchor["event"]))

    single = [(
        r["chrom"], r["pos"], r["cipos"],
        r["chrom"] in vector_set, r["qual"], r["vf"], r["sr"], r["homlen"],
        r["homseq"], r["n_beid"], r["imprecise"], r["sb"], r["event"], r["breakend"])
        for r in singles.values()]

    # Rank each table by GRIDSS QUAL (descending); qual "" sorts last.
    def _qual_key(row_qual):
        return row_qual if isinstance(row_qual, float) else float("-inf")

    integration.sort(key=lambda row: _qual_key(row[6]), reverse=True)
    vector_sv.sort(key=lambda row: _qual_key(row[5]), reverse=True)
    host_host.sort(key=lambda row: _qual_key(row[6]), reverse=True)
    single.sort(key=lambda row: _qual_key(row[4]), reverse=True)
    return integration, vector_sv, host_host, single


_INTEGRATION_HEADER = (
    "host_contig", "host_pos", "CIPOS",
    "vector_contig", "vector_pos", "breakend", "QUAL",
    "VF", "SR", "RP", "ASSR", "ASRP", "HOMLEN", "HOMSEQ",
    "n_beid", "IMPRECISE", "SB", "EVENT")
_INTEGRATION_HEADER_NLAIII = _INTEGRATION_HEADER + (
    "nlaiii_distance_host", "nlaiii_distance_vector", "nlaiii_category")

_HOST_HOST_HEADER = (
    "host_a", "pos_a", "CIPOS", "host_b", "pos_b", "breakend",
    "QUAL", "VF", "SR", "RP", "ASSR", "ASRP", "HOMLEN", "HOMSEQ",
    "n_beid", "IMPRECISE", "SB", "EVENT")
_HOST_HOST_HEADER_NLAIII = _HOST_HOST_HEADER + (
    "nlaiii_distance_1", "nlaiii_distance_2", "nlaiii_category")

_VECTOR_SV_HEADER = (
    "contig_a", "pos_a", "contig_b", "pos_b", "breakend", "QUAL",
    "VF", "SR", "RP", "ASSR", "ASRP", "HOMLEN", "HOMSEQ",
    "n_beid", "IMPRECISE", "SB", "EVENT")
_VECTOR_SV_HEADER_NLAIII = _VECTOR_SV_HEADER + (
    "nlaiii_distance_1", "nlaiii_distance_2", "nlaiii_category")
_SINGLE_HEADER = (
    "contig", "pos", "CIPOS", "is_vector", "QUAL",
    "VF", "SR", "HOMLEN", "HOMSEQ", "n_beid", "IMPRECISE", "SB", "EVENT", "breakend")

# Combined header for gridss_breakpoints.tsv: union of all four per-type
# headers. Each row carries a `breakpoint_type` prefix; columns absent for a
# given type are written as empty strings.
_BREAKPOINTS_BASE_HEADER = (
    "breakpoint_type",
    "host_contig", "host_pos", "CIPOS", "vector_contig", "vector_pos",
    "contig_a", "pos_a", "contig_b", "pos_b",
    "host_a", "host_b", "contig", "pos", "is_vector",
    "breakend", "QUAL", "VF", "SR", "RP", "ASSR", "ASRP",
    "HOMLEN", "HOMSEQ", "n_beid", "IMPRECISE", "SB", "EVENT",
)
_BREAKPOINTS_HEADER_NLAIII = _BREAKPOINTS_BASE_HEADER + (
    "nlaiii_distance_host", "nlaiii_distance_vector",
    "nlaiii_distance_1", "nlaiii_distance_2", "nlaiii_category",
)


def _count_breakpoints_by_type(path):
    """Return (type_counts, nlaiii_both_counts) dicts from a combined breakpoints TSV."""
    type_counts = {}
    nlaiii_both = {}
    with open(path) as fh:
        header = next(fh).rstrip("\n").split("\t")
        bt_idx = header.index("breakpoint_type") if "breakpoint_type" in header else None
        nl_idx = header.index("nlaiii_category") if "nlaiii_category" in header else None
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            bt = cols[bt_idx] if bt_idx is not None and bt_idx < len(cols) else ""
            type_counts[bt] = type_counts.get(bt, 0) + 1
            if nl_idx is not None and nl_idx < len(cols) and cols[nl_idx] == "both_ends_near_nlaiii":
                nlaiii_both[bt] = nlaiii_both.get(bt, 0) + 1
    return type_counts, nlaiii_both


def _write_breakpoints_tsv(path, integration, vector_sv, host_host, single,
                           integration_header, host_host_header, combined_header,
                           vector_sv_header=None):
    """Write all four breakpoint categories to one combined TSV."""
    if vector_sv_header is None:
        vector_sv_header = _VECTOR_SV_HEADER
    table_map = [
        ("integration_sites", integration, integration_header),
        ("vector_sv", vector_sv, vector_sv_header),
        ("host_host_bnd", host_host, host_host_header),
        ("single_breakend", single, _SINGLE_HEADER),
    ]
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(combined_header) + "\n")
        for btype, rows, hdr in table_map:
            for row in rows:
                d = {"breakpoint_type": btype}
                d.update(zip(hdr, row))
                fh.write("\t".join(str(d.get(col, "")) for col in combined_header) + "\n")
    os.replace(tmp, path)


def run(vcf_path, plasmid_contigs, out_prefix, mixed_fa=None):
    """Classify GRIDSS PASS breakends into a combined gridss_breakpoints.tsv.

    All four breakpoint categories (integration_sites, vector_sv, host_host_bnd,
    single_breakend) are written to one TSV with a ``breakpoint_type`` column.
    Skips writing when the combined file already exists. Returns
    ``{"outputs": {"breakpoints": path}, "counts": {...}}`` where counts include
    per-category row counts and, when NlaIII annotation is active, the
    both-ends-near-nlaiii counts for integration and host-host categories.

    When ``mixed_fa`` is provided, NlaIII (CATG) cut sites are computed and
    each host<->vector and host<->host junction row is annotated with
    breakend-to-nearest-cut distances and a category label. NlaIII counts are
    included in the returned counts dict for summary reporting. No GRIDSS calls
    are removed; the annotation is exploratory.
    """
    breakpoints_path = f"{out_prefix}.gridss_breakpoints.tsv"
    outputs = {"breakpoints": breakpoints_path}
    if is_nonempty(breakpoints_path):
        log.info(f"Breakpoints table already present, skipping: {breakpoints_path}")
        type_counts, nlaiii_both = _count_breakpoints_by_type(breakpoints_path)
        counts = {
            "integration_sites": type_counts.get("integration_sites", 0),
            "vector_sv": type_counts.get("vector_sv", 0),
            "host_host_bnd": type_counts.get("host_host_bnd", 0),
            "single_breakend": type_counts.get("single_breakend", 0),
            "n_host_vector_both_ends_nlaiii": nlaiii_both.get("integration_sites", 0),
            "n_host_host_both_ends_nlaiii": nlaiii_both.get("host_host_bnd", 0),
        }
        return {"outputs": outputs, "counts": counts}
    if not is_nonempty(vcf_path):
        raise FileNotFoundError(f"GRIDSS VCF missing or empty: {vcf_path}")

    vector_set = read_plasmid_contigs(plasmid_contigs)
    integration, vector_sv, host_host, single = _classify(
        _read_pass_records(vcf_path), vector_set)

    # NlaIII restriction-site annotation (exploratory; does not filter calls).
    # Both host<->vector and host<->host junctions are annotated. NlaIII counts
    # are included in the returned counts dict for summary reporting.
    cut_map = None
    if mixed_fa is not None:
        log.info(f"Building NlaIII (CATG) cut map from "
                 f"{os.path.basename(mixed_fa)} ...")
        cut_map = rsa.build_nlaiii_cut_map(mixed_fa)
        rsa.log_cut_map_summary(cut_map, vector_set)

    integration_header = _INTEGRATION_HEADER
    vector_sv_header = _VECTOR_SV_HEADER
    host_host_header = _HOST_HOST_HEADER
    combined_header = _BREAKPOINTS_BASE_HEADER
    if cut_map is not None:
        integration = rsa.annotate_junctions(integration, cut_map)
        vector_sv = rsa.annotate_vector_sv(vector_sv, cut_map)
        host_host = rsa.annotate_junctions(host_host, cut_map)
        integration_header = _INTEGRATION_HEADER_NLAIII
        vector_sv_header = _VECTOR_SV_HEADER_NLAIII
        host_host_header = _HOST_HOST_HEADER_NLAIII
        combined_header = _BREAKPOINTS_HEADER_NLAIII
        log.info(
            "NlaIII host-vector both/one/neither: "
            + "/".join(
                str(sum(1 for r in integration if r[-1] == c))
                for c in ("both_ends_near_nlaiii",
                          "one_end_near_nlaiii",
                          "neither_end_near_nlaiii")))
        log.info(
            "NlaIII host-host both/one/neither: "
            + "/".join(
                str(sum(1 for r in host_host if r[-1] == c))
                for c in ("both_ends_near_nlaiii",
                          "one_end_near_nlaiii",
                          "neither_end_near_nlaiii")))

    nlaiii_hv_both = (
        sum(1 for r in integration if r[-1] == "both_ends_near_nlaiii")
        if cut_map else 0)
    nlaiii_hh_both = (
        sum(1 for r in host_host if r[-1] == "both_ends_near_nlaiii")
        if cut_map else 0)

    _write_breakpoints_tsv(
        breakpoints_path, integration, vector_sv, host_host, single,
        integration_header, host_host_header, combined_header,
        vector_sv_header=vector_sv_header)

    counts = {
        "integration_sites": len(integration),
        "vector_sv": len(vector_sv),
        "host_host_bnd": len(host_host),
        "single_breakend": len(single),
        "n_host_vector_both_ends_nlaiii": nlaiii_hv_both,
        "n_host_host_both_ends_nlaiii": nlaiii_hh_both,
    }
    log.info(f"Classified {len(integration)} host<->vector, {len(vector_sv)} "
             f"vector<->vector, {len(host_host)} host<->host, {len(single)} "
             f"single-breakend junctions")
    log.info(f"  breakpoints -> {breakpoints_path}")
    return {"outputs": outputs, "counts": counts}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Classify GRIDSS PASS breakends into host<->vector "
                    "integration sites, vector<->vector SV candidates, "
                    "host<->host (manual-inspection), and single-breakend "
                    "tables.")
    parser.add_argument("-v", "--vcf", required=True,
                        help="GRIDSS breakend VCF (BGZF)")
    parser.add_argument("-p", "--plasmid_contigs", required=True,
                        help="Vector/plasmid contig names, one per line")
    parser.add_argument("-o", "--out_prefix", required=True,
                        help="Output path prefix for the four TSV tables")
    parser.add_argument("-r", "--mixed_fa", default=None,
                        help="Mixed reference FASTA (host+vector). When provided, "
                             "NlaIII (CATG) cut sites are annotated on each junction "
                             "and a restriction_site_gridss_summary.tsv is written.")
    args = parser.parse_args()
    run(args.vcf, args.plasmid_contigs, args.out_prefix, mixed_fa=args.mixed_fa)
