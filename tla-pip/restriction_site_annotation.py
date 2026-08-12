#!/usr/bin/env python3

"""NlaIII restriction-site annotation for TLA-seq GRIDSS junctions (v1).

TLA-seq libraries are built by NlaIII digestion (motif CATG, cuts 3' of the
motif on both strands — CATG is its own reverse complement, so one forward
scan covers both strands). This module annotates junction tables from
``gridss_sites.py`` with the distance from each breakend to the nearest
predicted NlaIII cut site in the mixed (host + vector) reference, and writes
an exploratory summary TSV. Proximity to a cut site is expected for any
correctly-formed restriction-ligation product (real or artifact); this
module only measures it, it does not itself call artifacts.

Coordinates are 1-based throughout (VCF POS convention, passed through
unmodified from ``gridss_sites.py``). A cut position is stored as the
1-based index of the motif's last base (0-based find offset + motif length).
"""

import bisect
from array import array

from pipeline_io import get_logger, write_tsv

log = get_logger()

NLAIII_MOTIF = "CATG"
NLAIII_PROXIMITY_BP = 2  # breakend counts as "near" a cut within this many bp
_MOTIF_LEN = len(NLAIII_MOTIF)

RESTRICTION_SUMMARY_HEADER = [
    "junction_type", "nlaiii_category", "count", "total", "percentage",
]


# ---------------------------------------------------------------------------
# Cut-site map
# ---------------------------------------------------------------------------

def build_nlaiii_cut_map(fasta_path):
    """Return {contig: array of sorted 1-based NlaIII cut positions}.

    Uses ``array('I', ...)`` rather than a plain list — at whole-genome
    scale (~1 CATG per 256bp by chance, so tens of millions of positions
    for a mammalian genome) this cuts memory roughly 5-7x versus a list of
    Python ints, and bisect works on it identically.
    """
    cut_map: dict[str, array] = {}
    with open(fasta_path) as fh:
        contig, chunks = None, []
        for line in fh:
            if line.startswith(">"):
                if contig is not None:
                    cut_map[contig] = _find_cuts("".join(chunks))
                contig, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip().upper())
        if contig is not None:
            cut_map[contig] = _find_cuts("".join(chunks))
    return cut_map


def _find_cuts(seq: str) -> array:
    """Sorted 1-based NlaIII cut positions in ``seq``.

    ``str.find`` in a loop, not regex: for a short literal motif this stays
    C-speed per call, which matters more than the Python-level loop overhead
    at genome scale. CATG can't overlap itself, so ``start`` jumps a full
    motif length past each hit rather than just +1.
    """
    positions = []
    start = 0
    while True:
        idx = seq.find(NLAIII_MOTIF, start)
        if idx == -1:
            break
        positions.append(idx + _MOTIF_LEN)
        start = idx + _MOTIF_LEN
    return array("I", positions)


def nearest_nlaiii_distance(cut_map, contig, pos):
    """Distance (bp, >=0) from 1-based ``pos`` on ``contig`` to nearest cut.

    None if ``contig`` isn't in ``cut_map`` or has no cut sites.
    """
    cuts = cut_map.get(contig)
    if not cuts:
        return None
    i = bisect.bisect_left(cuts, pos)
    candidates = []
    if i < len(cuts):
        candidates.append(abs(cuts[i] - pos))
    if i > 0:
        candidates.append(abs(cuts[i - 1] - pos))
    return min(candidates)


# ---------------------------------------------------------------------------
# Junction classification and annotation
# ---------------------------------------------------------------------------

def nlaiii_category(dist1, dist2, tolerance=NLAIII_PROXIMITY_BP):
    """Classify a junction by breakend proximity to NlaIII cut sites.

    ``dist1``/``dist2``: distance or None (no cuts on that contig; treated
    conservatively as "not near"). Returns 'both_ends_near_nlaiii',
    'one_end_near_nlaiii', or 'neither_end_near_nlaiii'.
    """
    near1 = dist1 is not None and dist1 <= tolerance
    near2 = dist2 is not None and dist2 <= tolerance
    if near1 and near2:
        return "both_ends_near_nlaiii"
    if near1 or near2:
        return "one_end_near_nlaiii"
    return "neither_end_near_nlaiii"


def annotate_vector_sv(vector_sv_rows, cut_map):
    """Annotate vector<->vector SV rows with NlaIII distances + category.

    Layout: (contig_a, pos_a, contig_b, pos_b, ...) — no CIPOS at index 2,
    so anchor 2 is at columns 2/3 (not 3/4 like integration/host_host).
    Appends (nlaiii_distance_1, nlaiii_distance_2, nlaiii_category).
    """
    annotated = []
    for row in vector_sv_rows:
        contig1, pos1 = row[0], row[1]
        contig2, pos2 = row[2], row[3]
        d1 = nearest_nlaiii_distance(cut_map, contig1, pos1)
        d2 = nearest_nlaiii_distance(cut_map, contig2, pos2)
        annotated.append(tuple(row) + (d1, d2, nlaiii_category(d1, d2)))
    return annotated


def annotate_junctions(junctions, cut_map):
    """Annotate two-sided junction rows with NlaIII distances + category.

    Works for both the integration (host<->vector) and host_host tables:
    both share the same row layout — anchor 1 at columns 0/1, anchor 2 at
    columns 3/4 (column 2 is CIPOS, unused here). Appends 3 columns:
    distance for anchor 1, distance for anchor 2, category. Returns a new
    list; does not mutate input.

    (Replaces the old separate ``annotate_integration_junctions`` /
    ``annotate_host_host_junctions`` — they had identical bodies since the
    row layout is shared, and ``annotate_integration_junctions`` took an
    unused ``vector_set`` argument.)
    """
    annotated = []
    for row in junctions:
        contig1, pos1 = row[0], row[1]
        contig2, pos2 = row[3], row[4]
        d1 = nearest_nlaiii_distance(cut_map, contig1, pos1)
        d2 = nearest_nlaiii_distance(cut_map, contig2, pos2)
        annotated.append(tuple(row) + (d1, d2, nlaiii_category(d1, d2)))
    return annotated


# ---------------------------------------------------------------------------
# Exploratory summary
# ---------------------------------------------------------------------------

def write_restriction_summary(path, integration_ann, host_host_ann):
    """Write category counts/percentages per junction type to a TSV."""
    cats = ("both_ends_near_nlaiii", "one_end_near_nlaiii", "neither_end_near_nlaiii")
    rows = []
    for jtype, data in (("host-vector", integration_ann), ("host-host", host_host_ann)):
        total = len(data)
        for cat in cats:
            count = sum(1 for row in data if row[-1] == cat)
            rows.append([jtype, cat, count, total, round(100.0 * count / total, 2) if total else 0.0])
    write_tsv(rows, RESTRICTION_SUMMARY_HEADER, path)
    log.info(f"NlaIII restriction-site summary -> {path}")


def log_cut_map_summary(cut_map, plasmid_set):
    """Log NlaIII cut-site counts in host vs. vector contigs."""
    host_cuts = sum(len(v) for k, v in cut_map.items() if k not in plasmid_set)
    vec_cuts = sum(len(v) for k, v in cut_map.items() if k in plasmid_set)
    log.info(f"NlaIII (CATG) cut sites: {host_cuts} in host contigs, {vec_cuts} in vector contigs")