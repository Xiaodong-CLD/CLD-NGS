#!/usr/bin/env python3

import argparse
import re

# A kept haplotype whose frequency over total merged reads is below this is
# marked flags=low_support (report with caution). Distinct from the reporting
# drop gate (--consensus_min_freq), which is pool-relative.
LOW_SUPPORT_FREQ = 0.01


def read_fasta(path):
    """Simple FASTA parser. Returns list of (header_without_>, sequence)."""
    records = []
    header = None
    seq_parts = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        records.append((header, "".join(seq_parts).upper()))

    return records


def get_size(header):
    """Extract VSEARCH size from header like seq2;size=18502;."""
    m = re.search(r"size=(\d+)", header)
    return int(m.group(1)) if m else "NA"


def clean_query_id(header):
    """
    Convert VSEARCH header:
      seq2;size=18502;
    into:
      seq2
    """
    first_token = header.split()[0]
    first_token = re.sub(r";size=\d+;?", "", first_token)
    first_token = first_token.rstrip(";")
    return first_token


def needleman_wunsch_edit_distance(ref, query, gap_open=2, gap_extend=1):
    """
    Global alignment using affine gap penalties (Gotoh):
      match    = 0
      mismatch = 1
      gap of length L = gap_open + L * gap_extend

    The gap-open penalty makes a single contiguous indel cheaper than several
    short indels threaded through free matched bases. In low-complexity tracts
    (e.g. the C/T-rich Fut8 66-78 window) this collapses what edit-distance
    scoring would scatter into one canonical deletion.

    Returns:
      ref_aln, query_aln, edit_distance
    where edit_distance is the plain mismatch + gap-base count (unchanged
    meaning), recomputed from the final alignment rather than the affine score.
    """

    n = len(ref)
    m = len(query)
    INF = float("inf")

    # M: ref[i-1] aligned to query[j-1]; X: ref base aligned to gap (deletion);
    # Y: query base aligned to gap (insertion).
    M = [[INF] * (m + 1) for _ in range(n + 1)]
    X = [[INF] * (m + 1) for _ in range(n + 1)]
    Y = [[INF] * (m + 1) for _ in range(n + 1)]
    tM = [[None] * (m + 1) for _ in range(n + 1)]
    tX = [[None] * (m + 1) for _ in range(n + 1)]
    tY = [[None] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0] = gap_open + i * gap_extend
        tX[i][0] = "X" if i > 1 else "M"
    for j in range(1, m + 1):
        Y[0][j] = gap_open + j * gap_extend
        tY[0][j] = "Y" if j > 1 else "M"

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            sub_cost = 0 if ref[i - 1] == query[j - 1] else 1

            best = min(
                (M[i - 1][j - 1], "M"),
                (X[i - 1][j - 1], "X"),
                (Y[i - 1][j - 1], "Y"),
                key=lambda t: t[0],
            )
            M[i][j] = best[0] + sub_cost
            tM[i][j] = best[1]

            bestX = min(
                (M[i - 1][j] + gap_open + gap_extend, "M"),
                (X[i - 1][j] + gap_extend, "X"),
                (Y[i - 1][j] + gap_open + gap_extend, "Y"),
                key=lambda t: t[0],
            )
            X[i][j] = bestX[0]
            tX[i][j] = bestX[1]

            bestY = min(
                (M[i][j - 1] + gap_open + gap_extend, "M"),
                (Y[i][j - 1] + gap_extend, "Y"),
                (X[i][j - 1] + gap_open + gap_extend, "X"),
                key=lambda t: t[0],
            )
            Y[i][j] = bestY[0]
            tY[i][j] = bestY[1]

    state = min(
        (M[n][m], "M"), (X[n][m], "X"), (Y[n][m], "Y"), key=lambda t: t[0]
    )[1]

    ref_aln = []
    query_aln = []
    i = n
    j = m

    while i > 0 or j > 0:
        if state == "M":
            ps = tM[i][j]
            ref_aln.append(ref[i - 1])
            query_aln.append(query[j - 1])
            i -= 1
            j -= 1
        elif state == "X":
            ps = tX[i][j]
            ref_aln.append(ref[i - 1])
            query_aln.append("-")
            i -= 1
        else:  # "Y"
            ps = tY[i][j]
            ref_aln.append("-")
            query_aln.append(query[j - 1])
            j -= 1
        state = ps

    ref_aln = "".join(reversed(ref_aln))
    query_aln = "".join(reversed(query_aln))

    edit_distance = sum(1 for r, q in zip(ref_aln, query_aln) if r != q)

    return ref_aln, query_aln, edit_distance


def call_variants_from_gapped_alignment(ref_aln, query_aln):
    """
    Convert gapped pairwise alignment into variant-style calls.

    Coordinates are 1-based reference coordinates.

    DEL:
      ref_start_1based = first deleted ref base
      ref_end_1based   = last deleted ref base

    INS:
      ref_start_1based = position before insertion
      ref_end_1based   = position before insertion
      If insertion is before first reference base, position is 0.

    SNV_or_MNV:
      one or more consecutive mismatches.
    """

    calls = []
    ref_pos = 0
    i = 0
    n = len(ref_aln)

    while i < n:
        r = ref_aln[i]
        q = query_aln[i]

        # Match
        if r != "-" and q != "-" and r == q:
            ref_pos += 1
            i += 1
            continue

        # Mismatch block: SNV or MNV
        if r != "-" and q != "-" and r != q:
            start_ref = ref_pos + 1
            ref_bases = []
            alt_bases = []

            while (
                i < n
                and ref_aln[i] != "-"
                and query_aln[i] != "-"
                and ref_aln[i] != query_aln[i]
            ):
                ref_bases.append(ref_aln[i])
                alt_bases.append(query_aln[i])
                ref_pos += 1
                i += 1

            calls.append({
                "type": "SNV_or_MNV",
                "ref_start_1based": start_ref,
                "ref_end_1based": ref_pos,
                "ref": "".join(ref_bases),
                "alt": "".join(alt_bases),
                "length_change": 0,
            })
            continue

        # Deletion relative to reference
        if r != "-" and q == "-":
            start_ref = ref_pos + 1
            deleted = []

            while i < n and ref_aln[i] != "-" and query_aln[i] == "-":
                deleted.append(ref_aln[i])
                ref_pos += 1
                i += 1

            calls.append({
                "type": "DEL",
                "ref_start_1based": start_ref,
                "ref_end_1based": ref_pos,
                "ref": "".join(deleted),
                "alt": "-",
                "length_change": -len(deleted),
            })
            continue

        # Insertion relative to reference
        if r == "-" and q != "-":
            insertion_anchor = ref_pos
            inserted = []

            while i < n and ref_aln[i] == "-" and query_aln[i] != "-":
                inserted.append(query_aln[i])
                i += 1

            calls.append({
                "type": "INS",
                "ref_start_1based": insertion_anchor,
                "ref_end_1based": insertion_anchor,
                "ref": "-",
                "alt": "".join(inserted),
                "length_change": len(inserted),
            })
            continue

        i += 1

    return calls


def is_edge_call(call, ref_len, edge_bp):
    """True if a variant lies within edge_bp of either reference end.

    These are masked from the consensus because amplicon ends carry primer-trim
    artifacts (e.g. the universal DEL 1 and DEL 129/130 in the Fut8 amplicon).
    """
    return call["ref_start_1based"] <= edge_bp or call["ref_end_1based"] > ref_len - edge_bp


def haplotype_signature(internal_calls):
    """Canonical string identifying a haplotype by its internal variants.

    Reads with the same signature are the same allele. WT (no internal
    variant) gets the literal signature "WT".
    """
    if not internal_calls:
        return "WT"
    parts = []
    for c in sorted(internal_calls, key=lambda c: c["ref_start_1based"]):
        allele = c["alt"] if c["type"] == "INS" else c["ref"]
        if c["type"] == "SNV_or_MNV":
            allele = f"{c['ref']}>{c['alt']}"
        parts.append(f"{c['ref_start_1based']}-{c['ref_end_1based']}:{c['type']}:{allele}")
    return "|".join(parts)


def reconstruct_consensus_ref_anchored(ref_seq, ref_aln, query_aln, edge_bp):
    """Build a reference-anchored consensus row for one haplotype's representative.

    The returned sequence has exactly len(ref_seq) characters -- one per
    reference position, 1:1 with the reference -- so every haplotype row and the
    reference line up column-for-column:

      - covered match / SNV                    -> the haplotype base
      - covered deletion                       -> '-' (ref base absent in haplotype)
      - within ``edge_bp`` of either ref end   -> the reference base (amplicon
        (and covered)                             ends carry primer-trim artifacts;
                                                  masked back to reference)
      - reference positions outside the         -> 'N' (the representative read did
        representative read's footprint            not span this region)

    Insertions (query bases with no reference position) cannot be shown in
    reference coordinates; they are returned separately as (anchor, seq) pairs,
    where ``anchor`` is the 1-based reference position the insertion follows
    (0 = before the first base).

    Returns (anchored_seq, cov_start_1based, cov_end_1based, insertions).
    """
    ref_len = len(ref_seq)
    out = ["N"] * ref_len

    # Covered column span = first..last alignment column carrying a query base.
    first_q = last_q = None
    for col, qb in enumerate(query_aln):
        if qb != "-":
            if first_q is None:
                first_q = col
            last_q = col
    if first_q is None:
        return "".join(out), 0, 0, []

    insertions = []          # (anchor_1based, inserted_seq)
    pending_anchor = None
    pending_seq = []
    cov_start = cov_end = None
    ref_pos = 0              # 0-based count of reference bases consumed so far

    for col in range(len(ref_aln)):
        rb = ref_aln[col]
        qb = query_aln[col]
        covered = first_q <= col <= last_q

        if rb != "-":
            p = ref_pos      # 0-based reference index
            pos1 = p + 1
            if not covered:
                out[p] = "N"
            else:
                if cov_start is None:
                    cov_start = pos1
                cov_end = pos1
                if pos1 <= edge_bp or pos1 > ref_len - edge_bp:
                    out[p] = ref_seq[p]         # edge mask -> reference base
                elif qb == "-":
                    out[p] = "-"                # deletion
                else:
                    out[p] = qb                 # match or SNV
            ref_pos += 1
        elif covered and qb != "-":
            # reference gap within the covered span = insertion
            if pending_anchor != ref_pos:
                if pending_seq:
                    insertions.append((pending_anchor, "".join(pending_seq)))
                pending_anchor = ref_pos
                pending_seq = []
            pending_seq.append(qb)

    if pending_seq:
        insertions.append((pending_anchor, "".join(pending_seq)))

    return "".join(out), (cov_start or 0), (cov_end or 0), insertions


def make_pretty_alignment(ref_aln, query_aln, width=80):
    lines = []

    for start in range(0, len(ref_aln), width):
        r = ref_aln[start:start + width]
        q = query_aln[start:start + width]

        mid = []
        for rb, qb in zip(r, q):
            if rb == qb:
                mid.append("|")
            elif rb == "-" or qb == "-":
                mid.append(" ")
            else:
                mid.append("*")

        lines.append(f"REF   {r}")
        lines.append(f"      {''.join(mid)}")
        lines.append(f"QUERY {q}")
        lines.append("")

    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", required=True, help="Reference amplicon FASTA")
    parser.add_argument("--query", required=True, help="Top VSEARCH sequences FASTA")
    parser.add_argument("--out_prefix", required=True, help="Output prefix for the consensus FASTA deliverable")
    parser.add_argument(
        "--diag_prefix", default=None,
        help="Output prefix for the diagnostic files (alignments.txt, "
             "variants.tsv, summary.tsv). Defaults to --out_prefix."
    )
    parser.add_argument(
        "--include_ref_match",
        action="store_true",
        help="If set, output REF_MATCH rows for sequences identical to reference."
    )
    parser.add_argument(
        "--consensus_edge_bp", type=int, default=5,
        help="Mask variants within this many bp of either reference end when "
             "building consensus haplotypes (primer-trim artifacts). Default 5."
    )
    parser.add_argument(
        "--consensus_min_size", type=int, default=50,
        help="Minimum summed read support for a haplotype to be reported. Default 50."
    )
    parser.add_argument(
        "--consensus_min_freq", type=float, default=0.01,
        help="Reporting drop gate: a haplotype whose pool-relative frequency "
             "(read_support / summed top-N support) is below this is not "
             "reported. Default 0.01 (1%%); set 0 to report all."
    )
    parser.add_argument(
        "--total_reads", type=int, default=0,
        help="Total merged reads for the sample, used as the displayed "
             "frequency denominator and the low_support flag basis. When 0 "
             "(unset), falls back to the summed read support of the input "
             "(top-N) sequences."
    )
    args = parser.parse_args(argv)

    ref_records = read_fasta(args.ref)
    if len(ref_records) != 1:
        raise ValueError("Reference FASTA should contain exactly one sequence.")

    ref_id, ref_seq = ref_records[0]

    query_records = read_fasta(args.query)

    variant_rows = []
    summary_rows = []
    # Per-haplotype accumulation keyed by edge-masked variant signature.
    haplotypes = {}
    total_size = 0

    diag_prefix = args.diag_prefix or args.out_prefix
    alignment_txt = diag_prefix + ".alignments.txt"
    variants_tsv = diag_prefix + ".variants.tsv"
    summary_tsv = diag_prefix + ".summary.tsv"
    consensus_fa = args.out_prefix + ".consensus.fa"

    with open(alignment_txt, "w") as aln_out:
        for header, query_seq in query_records:
            query_id = clean_query_id(header)
            size = get_size(header)

            ref_aln, query_aln, edit_distance = needleman_wunsch_edit_distance(
                ref_seq,
                query_seq
            )

            calls = call_variants_from_gapped_alignment(ref_aln, query_aln)

            # Group into a haplotype by its internal (edge-masked) variants.
            internal_calls = [
                c for c in calls
                if not is_edge_call(c, len(ref_seq), args.consensus_edge_bp)
            ]
            sig = haplotype_signature(internal_calls)
            size_int = size if isinstance(size, int) else 0
            total_size += size_int
            hap = haplotypes.get(sig)
            if hap is None:
                haplotypes[sig] = {
                    "size": size_int,
                    "n_members": 1,
                    "internal_calls": internal_calls,
                    "top_size": size_int,
                    "rep_query_id": query_id,
                    "rep_ref_aln": ref_aln,
                    "rep_query_aln": query_aln,
                }
            else:
                hap["size"] += size_int
                hap["n_members"] += 1
                if size_int > hap["top_size"]:
                    hap["top_size"] = size_int
                    hap["internal_calls"] = internal_calls
                    hap["rep_query_id"] = query_id
                    hap["rep_ref_aln"] = ref_aln
                    hap["rep_query_aln"] = query_aln

            n_mismatch = sum(1 for r, q in zip(ref_aln, query_aln)
                             if r != "-" and q != "-" and r != q)
            n_ins = sum(1 for r, q in zip(ref_aln, query_aln)
                        if r == "-" and q != "-")
            n_del = sum(1 for r, q in zip(ref_aln, query_aln)
                        if r != "-" and q == "-")

            summary_rows.append([
                query_id,
                size,
                len(query_seq),
                edit_distance,
                n_mismatch,
                n_ins,
                n_del,
                len(calls),
            ])

            if calls:
                for variant_id, call in enumerate(calls, start=1):
                    variant_rows.append([
                        query_id,
                        size,
                        len(query_seq),
                        variant_id,
                        call["type"],
                        call["ref_start_1based"],
                        call["ref_end_1based"],
                        call["ref"],
                        call["alt"],
                        call["length_change"],
                    ])
            elif args.include_ref_match:
                variant_rows.append([
                    query_id,
                    size,
                    len(query_seq),
                    0,
                    "REF_MATCH",
                    ".",
                    ".",
                    ".",
                    ".",
                    0,
                ])

            aln_out.write(f">{header}\n")
            aln_out.write(f"query_id={query_id}\n")
            aln_out.write(f"size={size}\n")
            aln_out.write(f"query_length={len(query_seq)}\n")
            aln_out.write(f"edit_distance={edit_distance}\n")
            aln_out.write(f"num_variant_events={len(calls)}\n\n")
            aln_out.write(make_pretty_alignment(ref_aln, query_aln))
            aln_out.write("\n" + "=" * 100 + "\n\n")

    with open(variants_tsv, "w") as out:
        out.write(
            "query_id\tread_support\tquery_length\tvariant_id\ttype\t"
            "ref_start_1based\tref_end_1based\tref\talt\tlength_change\n"
        )
        for row in variant_rows:
            out.write("\t".join(map(str, row)) + "\n")

    with open(summary_tsv, "w") as out:
        out.write(
            "query_id\tread_support\tquery_length\tedit_distance\t"
            "mismatches\tinserted_bases\tdeleted_bases\tn_variant_events\n"
        )
        for row in summary_rows:
            out.write("\t".join(map(str, row)) + "\n")

    # Consensus haplotypes: reference-anchored, frequency-sorted FASTA. The
    # reference is written first as a full-length baseline the haplotype rows
    # line up against. A haplotype is dropped below the absolute support floor
    # (--consensus_min_size) or the pool-relative frequency gate
    # (--consensus_min_freq, denominator = summed top-N support). The DISPLAYED
    # freq uses total merged reads (--total_reads, falling back to the pool);
    # a kept haplotype whose total-read freq is below LOW_SUPPORT_FREQ is
    # reported with flags=low_support rather than dropped.
    total_reads = args.total_reads if args.total_reads > 0 else total_size
    ranked = sorted(haplotypes.items(), key=lambda kv: kv[1]["size"], reverse=True)
    n_written = 0
    with open(consensus_fa, "w") as out:
        out.write(f">reference id={ref_id} length={len(ref_seq)}\n")
        out.write(ref_seq + "\n")
        for rank, (sig, hap) in enumerate(ranked, start=1):
            pool_freq = hap["size"] / total_size if total_size else 0.0
            if hap["size"] < args.consensus_min_size or pool_freq < args.consensus_min_freq:
                continue
            freq = hap["size"] / total_reads if total_reads else 0.0
            anchored_seq, cov_start, cov_end, insertions = reconstruct_consensus_ref_anchored(
                ref_seq, hap["rep_ref_aln"], hap["rep_query_aln"], args.consensus_edge_bp
            )
            ins_str = ",".join(f"after{a}:{s}" for a, s in insertions) if insertions else "none"
            flags = "low_support" if freq < LOW_SUPPORT_FREQ else "none"
            out.write(
                f">hap{rank}; read_support={hap['size']}; freq={freq:.6f}; "
                f"covered={cov_start}-{cov_end}; "
                f"n_internal_variants={len(hap['internal_calls'])}; "
                f"insertions={ins_str}; sig={sig}; flags={flags}\n"
            )
            out.write(anchored_seq + "\n")
            n_written += 1

    print(f"Wrote: {variants_tsv}")
    print(f"Wrote: {summary_tsv}")
    print(f"Wrote: {alignment_txt}")
    print(f"Wrote: {consensus_fa} ({n_written} haplotype(s))")


if __name__ == "__main__":
    main()