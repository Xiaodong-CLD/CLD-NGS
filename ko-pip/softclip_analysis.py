#!/usr/bin/env python3
"""Extract, cluster, and BLAST soft-clipped sequences in one pass.

Pipeline:
  1. Walk BAM and extract soft-clipped segments (5' and 3')
  2. Cluster near-identical clips with starcode (per side)
  3. BLAST top-N cluster representatives against 1+ local plasmid subjects;
     for each cluster keep the best hit per subject (top N = #subjects)
  4. Write cluster table (with blast_result + flag) and the BLAST report

Per top cluster, blast_result is 'hit' (pident >= 90% and coverage >= 80% of
the clip), 'partial' (aligns but below that), or 'no hit'. The flag column is
'construct match' (hit), 'insert like' (no hit), or blank, and is only set for
clusters with support_pct >= 1% of total analyzed reads AND clip length >= 30
bp; partial hits never flag. support_pct = cluster reads / total analyzed
reads (all reads passing the flag filter in the input BAM).

Intended for: breakpoint analysis on the full sample BAM (soft-clip junctions
left by large insertions / construct integrations).
"""

from __future__ import annotations

import argparse
import hashlib
import pathlib
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from typing import Dict, Iterator, List, Tuple

import pysam

LOG_PREFIX = "[Takeda Pipeline Message]"
EXCLUDE_FLAGS = 4 | 256 | 1024 | 2048  # 3332: unmapped/secondary/dup/supp
CIGAR_S = 4
BLAST_DB_EXTS = (".nhr", ".nin", ".nsq")

# Insert-flag thresholds (see classify_clusters).
FULL_HIT_MIN_PIDENT = 90.0   # percent identity for a "hit" (construct present)
FULL_HIT_MIN_COV = 80.0      # percent of the clip that must align for a "hit"
MIN_SUPPORT_FRAC = 0.01      # cluster reads / total analyzed reads (full BAM)
MIN_FLAG_CLIP_LEN = 30       # minimum representative clip length (bp) to flag


def log(msg: str) -> None:
    print(f"{LOG_PREFIX} {msg}", flush=True)


def count_total_reads(bam_path: pathlib.Path) -> int:
    """Count reads passing the EXCLUDE_FLAGS filter (the flag denominator)."""
    n = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.flag & EXCLUDE_FLAGS:
                continue
            n += 1
    return n



# ============================================================================
# EXTRACTION PHASE
# ============================================================================

def iter_softclips(
    bam_path: pathlib.Path,
    min_len: int,
) -> Iterator[Tuple[str, str, str, int]]:
    """Yield (clip_seq, side, chrom, ref_anchor_1based) per soft clip."""
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.flag & EXCLUDE_FLAGS:
                continue
            cigar = read.cigartuples
            seq = read.query_sequence
            if not cigar or not seq:
                continue
            chrom = read.reference_name or "*"
            # 5' clip
            op, length = cigar[0]
            if op == CIGAR_S and length >= min_len:
                clip = seq[:length]
                anchor = (read.reference_start + 1) if read.reference_start is not None else 0
                yield clip, "5p", chrom, anchor
            # 3' clip
            op, length = cigar[-1]
            if op == CIGAR_S and length >= min_len:
                clip = seq[-length:]
                anchor = read.reference_end if read.reference_end is not None else 0  # 1-based last aligned
                yield clip, "3p", chrom, anchor


def _run_starcode(
    seq_counts: Counter,
    edit_dist: int,
    algorithm: str,
    threads: int,
    workdir: pathlib.Path,
    side_tag: str,
) -> List[Tuple[str, int, List[str]]]:
    """Run starcode clustering and return (centroid, total_count, members) tuples."""
    if not seq_counts:
        return []

    starcode = shutil.which("starcode")
    if starcode is None:
        sys.exit(
            f"{LOG_PREFIX} ERROR: starcode binary not on PATH. "
            "Install with `conda install -n ko_pip_dev -c bioconda starcode`."
        )

    in_path = workdir / f"{side_tag}.counts.tsv"
    out_path = workdir / f"{side_tag}.clusters.tsv"
    with in_path.open("w") as fh:
        for seq, n in seq_counts.items():
            fh.write(f"{seq}\t{n}\n")

    cmd = [
        starcode,
        "-i", str(in_path),
        "-o", str(out_path),
        "-d", str(edit_dist),
        "-t", str(threads),
        "--print-clusters",
        "--quiet",
    ]
    if algorithm == "sphere":
        cmd.append("--sphere")
    elif algorithm == "connected":
        cmd.append("--connected-comp")

    log(f"  starcode {side_tag}: {len(seq_counts)} unique seqs, d={edit_dist}, algo={algorithm}")
    subprocess.run(cmd, check=True)

    clusters: List[Tuple[str, int, List[str]]] = []
    with out_path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            centroid, count_str, members_str = parts[0], parts[1], parts[2]
            members = members_str.split(",") if members_str else [centroid]
            clusters.append((centroid, int(count_str), members))
    return clusters


def cluster_with_starcode(
    records: List[Tuple[str, str, str, int]],
    edit_dist: int,
    algorithm: str,
    threads: int,
) -> List[Dict]:
    """Cluster soft-clip records via starcode, separately per side."""
    seq_counts_by_side: Dict[str, Counter] = {"5p": Counter(), "3p": Counter()}
    anchors_by_seq: Dict[Tuple[str, str], Counter] = {}

    for clip, side, chrom, anchor in records:
        seq_counts_by_side[side][clip] += 1
        key = (side, clip)
        anchors_by_seq.setdefault(key, Counter())[(chrom, anchor)] += 1

    out: List[Dict] = []
    with tempfile.TemporaryDirectory(prefix="softclip_starcode_") as td:
        workdir = pathlib.Path(td)
        for side in ("5p", "3p"):
            sc = seq_counts_by_side[side]
            if not sc:
                continue
            for centroid, total, members in _run_starcode(
                sc, edit_dist, algorithm, threads, workdir, side_tag=side
            ):
                cluster_anchors: Counter = Counter()
                for m in members:
                    cluster_anchors.update(anchors_by_seq.get((side, m), {}))
                (top_chrom, top_anchor), loc_n = cluster_anchors.most_common(1)[0]
                out.append({
                    "side": side,
                    "reads": total,
                    "members": len(members),
                    "rep_seq": centroid,
                    "rep_reads": sc.get(centroid, 0),
                    "top_chrom": top_chrom,
                    "top_anchor": top_anchor,
                    "locus_reads": loc_n,
                })

    out.sort(key=lambda c: -c["reads"])
    return out


def _write_top_clusters_tsv(
    clusters: List[Dict],
    out_dir: pathlib.Path,
    top_n: int,
    total_reads: int,
    edit_dist: int = 8,
    algorithm: str = "sphere",
    min_len: int = 10,
    loose_criteria_used: bool = False,
    classifications: Dict[int, dict] | None = None,
) -> pathlib.Path:
    """Write top_clusters.tsv with metadata header and top-N clusters.

    Metadata is written as comment lines (# prefix) at the top. ``support_pct``
    is the cluster's read count over ``total_reads`` (all analyzed reads in the
    BAM, post flag-filter) -- an event-fraction proxy. When ``classifications``
    is provided, the per-cluster ``blast_result`` and ``flag`` columns are
    appended (else '.' / '').
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = out_dir / "top_clusters.tsv"
    denom = total_reads if total_reads else 1
    classifications = classifications or {}

    with tsv_path.open("w") as fh:
        # Write metadata as comment lines
        fh.write("# Softclip Clustering Settings\n")
        fh.write(f"# edit_dist: {edit_dist}\n")
        fh.write(f"# algorithm: {algorithm}\n")
        fh.write(f"# min_len: {min_len}\n")
        fh.write(f"# loose_criteria_used: {loose_criteria_used}\n")
        fh.write(f"# total_analyzed_reads: {total_reads}\n")
        fh.write(f"# min_support_pct: {MIN_SUPPORT_FRAC * 100:.3f}\n")
        # Write header (support_pct = reads / total_analyzed_reads)
        fh.write("rank\treads\tsupport_pct\tside\tlength\ttop_locus\trep_sequence\tblast_result\tflag\n")
        # Write top-N clusters
        for rank, c in enumerate(clusters[:top_n], start=1):
            support_pct = 100.0 * c["reads"] / denom
            cls = classifications.get(rank, {})
            blast_result = cls.get("blast_result", ".")
            flag = cls.get("flag", "")
            fh.write(
                f"{rank}\t{c['reads']}\t{support_pct:.3f}\t{c['side']}\t{len(c['rep_seq'])}\t"
                f"{c['top_chrom']}:{c['top_anchor']}\t{c['rep_seq']}\t{blast_result}\t{flag}\n"
            )
    log(f"Wrote {tsv_path} ({min(top_n, len(clusters))} clusters, {len(clusters)} total)")
    return tsv_path



def write_cluster_table_and_fasta(
    clusters: List[Dict],
    out_dir: pathlib.Path,
    top_n: int,
    write_tsv: bool = False,
) -> Tuple[pathlib.Path, pathlib.Path]:
    """Write top_softclip_sequence.fa (and optionally softclip_clusters.tsv).
    
    Returns (clusters_tsv or None, fasta_path).
    Clusters info is still needed for BLAST summary, so returned as internal dict.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = out_dir / "softclip_clusters.tsv"
    fasta_path = out_dir / "top_softclip_sequence.fa"
    total = sum(c["reads"] for c in clusters)
    denom = total if total else 1

    # Optionally write TSV (can be large)
    if write_tsv:
        with tsv_path.open("w") as fh:
            fh.write(
                "rank\treads\tfreq_pct\tmembers\trep_reads\tside\tlength\t"
                "top_locus\tlocus_reads\trep_sequence\n"
            )
            for rank, c in enumerate(clusters, start=1):
                freq = 100.0 * c["reads"] / denom
                fh.write(
                    f"{rank}\t{c['reads']}\t{freq:.3f}\t{c['members']}\t{c['rep_reads']}\t"
                    f"{c['side']}\t{len(c['rep_seq'])}\t"
                    f"{c['top_chrom']}:{c['top_anchor']}\t{c['locus_reads']}\t{c['rep_seq']}\n"
                )
        log(f"Wrote {tsv_path}")

    with fasta_path.open("w") as fh:
        for rank, c in enumerate(clusters[: top_n], start=1):
            freq = 100.0 * c["reads"] / denom
            header = (
                f">cluster{rank} side={c['side']} reads={c['reads']} freq={freq:.2f}% "
                f"members={c['members']} rep_reads={c['rep_reads']} "
                f"len={len(c['rep_seq'])} "
                f"top_locus={c['top_chrom']}:{c['top_anchor']} locus_reads={c['locus_reads']}"
            )
            fh.write(f"{header}\n{c['rep_seq']}\n")

    log(f"Wrote {fasta_path} (top {min(top_n, len(clusters))})")
    return tsv_path, fasta_path


# ============================================================================
# BLAST PHASE
# ============================================================================

def _have_blast_db(db_prefix: pathlib.Path) -> bool:
    """True if all required BLAST DB files exist."""
    return all(pathlib.Path(str(db_prefix) + ext).exists() for ext in BLAST_DB_EXTS)


def _db_is_stale(db_prefix: pathlib.Path, sources: List[pathlib.Path]) -> bool:
    """True if any source FASTA is newer than DB files."""
    db_mtime = min(
        pathlib.Path(str(db_prefix) + ext).stat().st_mtime for ext in BLAST_DB_EXTS
    )
    return any(s.stat().st_mtime > db_mtime for s in sources)


def _combined_db_prefix(db_dir: pathlib.Path, subjects: List[pathlib.Path]) -> pathlib.Path:
    """Deterministic DB prefix derived from sorted absolute subject paths."""
    key = "\n".join(sorted(str(s.resolve()) for s in subjects)).encode()
    digest = hashlib.sha1(key).hexdigest()[:10]
    return db_dir / f"subjects_{digest}"


def prepare_blast_db(
    subjects: List[pathlib.Path],
    db_dir: pathlib.Path,
) -> pathlib.Path:
    """Return a ready-to-use BLAST DB prefix; build it if needed."""
    if len(subjects) == 1 and _have_blast_db(subjects[0].with_suffix("")):
        prefix = subjects[0].with_suffix("")
        log(f"Reusing existing BLAST DB at {prefix}.*")
        return prefix
    if len(subjects) == 1 and _have_blast_db(subjects[0]):
        log(f"Reusing existing BLAST DB at {subjects[0]}.*")
        return subjects[0]

    db_prefix = _combined_db_prefix(db_dir, subjects)
    db_prefix.parent.mkdir(parents=True, exist_ok=True)

    if _have_blast_db(db_prefix) and not _db_is_stale(db_prefix, subjects):
        log(f"Reusing cached BLAST DB at {db_prefix}.*")
        return db_prefix

    combined = pathlib.Path(str(db_prefix) + ".fa")
    log(f"Building BLAST DB at {db_prefix}.* from {len(subjects)} subject FASTA(s)")
    with combined.open("w") as out_fh:
        for s in subjects:
            with s.open() as in_fh:
                shutil.copyfileobj(in_fh, out_fh)
                out_fh.write("\n")

    makeblastdb = shutil.which("makeblastdb")
    if makeblastdb is None:
        sys.exit(
            f"{LOG_PREFIX} ERROR: makeblastdb not on PATH. "
            "Install with `conda install -n ko_pip_dev -c bioconda blast`."
        )
    subprocess.run(
        [makeblastdb, "-in", str(combined), "-dbtype", "nucl",
         "-out", str(db_prefix), "-title", db_prefix.name],
        check=True,
    )
    return db_prefix


def _parse_subjects(arg: str) -> List[pathlib.Path]:
    """Parse comma-separated subject paths."""
    paths = [pathlib.Path(p.strip()) for p in arg.split(",") if p.strip()]
    for p in paths:
        if not p.exists():
            sys.exit(f"{LOG_PREFIX} ERROR: subject FASTA not found: {p}")
    return paths


BLAST_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "sstrand",
]
BLAST_OUTFMT = "6 " + " ".join(BLAST_FIELDS)


def run_blastn(
    query_fasta: pathlib.Path,
    db_prefix: pathlib.Path,
    threads: int,
) -> List[dict]:
    """Run local blastn and return parsed hits.
    
    Uses BLAST's standard defaults: e-value=1e-5, word_size=11.
    """
    blastn = shutil.which("blastn")
    if blastn is None:
        sys.exit(
            f"{LOG_PREFIX} ERROR: blastn not on PATH. "
            "Install with `conda install -n ko_pip_dev -c bioconda blast`."
        )
    log(f"blastn {query_fasta.name} vs {pathlib.Path(str(db_prefix)).name} "
        f"(threads={threads})")
    proc = subprocess.run(
        [
            blastn,
            "-query", str(query_fasta),
            "-db", str(db_prefix),
            "-outfmt", BLAST_OUTFMT,
            "-num_threads", str(threads),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    int_cols = {"length", "mismatch", "gapopen", "qstart", "qend",
                "sstart", "send", "qlen", "slen"}
    float_cols = {"pident", "evalue", "bitscore"}
    rows: List[dict] = []
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != len(BLAST_FIELDS):
            continue
        row = dict(zip(BLAST_FIELDS, parts))
        for k in int_cols:
            row[k] = int(row[k])
        for k in float_cols:
            row[k] = float(row[k])
        rows.append(row)
    return rows


def _read_query_ids_and_lengths(fasta: pathlib.Path) -> Tuple[List[str], dict]:
    """Return (ordered query IDs, {id: length}) from FASTA."""
    order: List[str] = []
    lengths: dict = {}
    current_id: str | None = None
    current_len = 0
    with fasta.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                current_id = line[1:].split()[0]
                order.append(current_id)
                current_len = 0
            else:
                current_len += len(line)
        if current_id is not None:
            lengths[current_id] = current_len
    return order, lengths


def classify_clusters(
    clusters: List[Dict],
    hits: List[dict],
    query_lengths: dict,
    n_subjects: int,
    total_analyzed_reads: int,
) -> Dict[int, dict]:
    """Classify each top cluster's BLAST outcome and insert flag.

    For every cluster (rank = 1-based position in the reads-sorted list) keep
    the best hit per subject sequence, then the top ``n_subjects`` of those.
    blast_result: 'hit' (any kept hit with pident >= FULL_HIT_MIN_PIDENT and
    coverage >= FULL_HIT_MIN_COV), else 'partial' (some alignment), else
    'no hit'; '.' when the cluster was not BLASTed.

    flag (gated by support >= MIN_SUPPORT_FRAC of total analyzed reads AND
    representative clip length >= MIN_FLAG_CLIP_LEN): 'construct match' when
    blast_result == 'hit', 'insert like' when blast_result == 'no hit'; blank
    otherwise (partial never flags).
    """
    from collections import defaultdict
    by_query: Dict[str, List[dict]] = defaultdict(list)
    for h in hits:
        by_query[h["qseqid"]].append(h)

    denom = total_analyzed_reads if total_analyzed_reads else 1
    out: Dict[int, dict] = {}
    for idx, c in enumerate(clusters, start=1):
        qid = f"cluster{idx}"
        qlen = query_lengths.get(qid) or len(c["rep_seq"])

        best_per_subj: Dict[str, dict] = {}
        for h in by_query.get(qid, []):
            s = h["sseqid"]
            if s not in best_per_subj or h["bitscore"] > best_per_subj[s]["bitscore"]:
                best_per_subj[s] = h
        top_hits = sorted(best_per_subj.values(), key=lambda r: -r["bitscore"])[:max(1, n_subjects)]

        if qid not in query_lengths:
            blast_result = "."
        elif not top_hits:
            blast_result = "no hit"
        else:
            is_full = any(
                h["pident"] >= FULL_HIT_MIN_PIDENT
                and (100.0 * h["length"] / qlen if qlen else 0.0) >= FULL_HIT_MIN_COV
                for h in top_hits
            )
            blast_result = "hit" if is_full else "partial"

        support_frac = c["reads"] / denom
        gated = support_frac >= MIN_SUPPORT_FRAC and len(c["rep_seq"]) >= MIN_FLAG_CLIP_LEN
        if gated and blast_result == "hit":
            flag = "construct match"
        elif gated and blast_result == "no hit":
            flag = "insert like"
        else:
            flag = ""

        out[idx] = {"blast_result": blast_result, "flag": flag,
                    "top_hits": top_hits, "qlen": qlen}
    return out


def write_blast_report(
    clusters: List[Dict],
    classifications: Dict[int, dict],
    out_dir: pathlib.Path,
    total_reads: int,
) -> pathlib.Path:
    """Write blast_hits_full.tsv: top-N (N = #subjects) hits per cluster.

    Always written (one row even for 'no hit' clusters) so a reviewer can see
    every top cluster's BLAST outcome. ``support_pct`` is reads / total_reads
    (all analyzed reads, post flag-filter).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    blast_tsv = out_dir / "blast_hits_full.tsv"
    header = [
        "cluster_id", "n_reads", "support_pct", "repseq_len",
        "blast_result", "blast_hit_sseqid", "pident", "evalue",
        "aln_len", "cov_pct",
    ]
    denom = total_reads if total_reads else 1

    with blast_tsv.open("w") as fh:
        fh.write("\t".join(header) + "\n")
        for idx, c in enumerate(clusters, start=1):
            cls = classifications.get(idx)
            if cls is None:
                continue
            qid = f"cluster{idx}"
            support_pct = 100.0 * c["reads"] / denom
            qlen = cls["qlen"] or len(c["rep_seq"])
            repseq_len = len(c["rep_seq"])
            if cls["top_hits"]:
                for h in cls["top_hits"]:
                    cov = 100.0 * h["length"] / qlen if qlen else 0.0
                    fh.write("\t".join([
                        qid, str(c["reads"]), f"{support_pct:.3f}", str(repseq_len),
                        cls["blast_result"], h["sseqid"], f"{h['pident']:.2f}",
                        f"{h['evalue']:.2e}", str(h["length"]), f"{cov:.1f}",
                    ]) + "\n")
            else:
                fh.write("\t".join([
                    qid, str(c["reads"]), f"{support_pct:.3f}", str(repseq_len),
                    cls["blast_result"], ".", ".", ".", ".", ".",
                ]) + "\n")
    log(f"Wrote {blast_tsv} ({len(classifications)} clusters classified)")
    return blast_tsv



def _top_n_fasta(src: pathlib.Path, n: int, dest: pathlib.Path) -> Tuple[pathlib.Path, int]:
    """Write the first n FASTA records of src to dest. Returns (dest, n_written)."""
    if n <= 0:
        return src, -1
    written = 0
    keep = False
    with src.open() as in_fh, dest.open("w") as out_fh:
        for line in in_fh:
            if line.startswith(">"):
                if written >= n:
                    break
                keep = True
                written += 1
            if keep:
                out_fh.write(line)
    return dest, written


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        description="Extract, cluster, and BLAST soft-clipped sequences."
    )
    ap.add_argument("-b", "--bam", required=True, type=pathlib.Path,
                    help="Input BAM (uninformative or full sample).")
    ap.add_argument("-o", "--out-dir", required=True, type=pathlib.Path,
                    help="Output directory for cluster and BLAST results.")
    
    # Extraction knobs
    ap.add_argument("--min-len", type=int, default=10,
                    help="Minimum soft-clip length in bp (default: 10).")
    ap.add_argument("--top", type=int, default=20,
                    help="Number of top clusters to represent in FASTA (default: 20).")
    ap.add_argument("--edit-dist", type=int, default=8,
                    help="Starcode Levenshtein distance (default: 8).")
    ap.add_argument("--algorithm", choices=("message-passing", "sphere", "connected"),
                    default="sphere", help="Starcode algorithm (default: sphere).")
    ap.add_argument("--threads", type=int, default=4,
                    help="Threads for starcode and blastn (default: 4).")
    
    # BLAST knobs
    ap.add_argument("--blast-subject", default=None,
                    help="Local subject FASTA (comma-separated list) or omit to skip BLAST.")
    ap.add_argument("--blast-db-dir", type=pathlib.Path, default=None,
                    help="Directory for BLAST DB cache (default: <out-dir>/blast_db).")
    ap.add_argument("--blast-top-n", type=int, default=20,
                    help="Number of clusters to BLAST (default: 20; 0=all).")
    args = ap.parse_args(argv)

    if not args.bam.exists():
        sys.exit(f"{LOG_PREFIX} ERROR: BAM not found: {args.bam}")

    # ---- STAGE 1: Extract and cluster ----
    log(f"Scanning soft clips in {args.bam} "
        f"(min_len={args.min_len}, edit_dist={args.edit_dist}, algo={args.algorithm})")
    records: List[Tuple[str, str, str, int]] = list(iter_softclips(args.bam, args.min_len))
    total = len(records)
    log(f"Collected {total} soft-clip segments; clustering with starcode…")

    clusters = cluster_with_starcode(
        records, edit_dist=args.edit_dist, algorithm=args.algorithm, threads=args.threads,
    )
    log(f"Formed {len(clusters)} clusters")
    
    # Retry with looser criteria if no clusters found
    loose_used = False
    if len(clusters) == 0 and len(records) > 0:
        log(f"No clusters with edit_dist={args.edit_dist}; retrying with edit_dist=12...")
        clusters = cluster_with_starcode(
            records, edit_dist=12, algorithm=args.algorithm, threads=args.threads,
        )
        log(f"Formed {len(clusters)} clusters with loose criteria")
        loose_used = True
    
    if len(clusters) == 0:
        log("No clusters found even with loose criteria; skipping output files")
        return

    total_analyzed_reads = count_total_reads(args.bam)
    log(f"Total reads (post-filter) for support denominator = {total_analyzed_reads}")

    # Generate FASTA of cluster representatives (BLAST queries).
    clusters_tsv, fasta_path = write_cluster_table_and_fasta(clusters, args.out_dir, args.top, write_tsv=False)

    # ---- STAGE 2: BLAST + insert classification ----
    subjects = _parse_subjects(args.blast_subject or "")
    classifications: Dict[int, dict] = {}
    if subjects:
        blast_db_dir = args.blast_db_dir or (args.out_dir / "blast_db")
        db_prefix = prepare_blast_db(subjects, blast_db_dir)

        # Restrict BLAST to the top-N cluster representatives.
        work_query, n_kept = _top_n_fasta(
            fasta_path, args.blast_top_n, args.out_dir / "blast_query_top.fasta"
        )
        if n_kept >= 0:
            log(f"Restricting BLAST to top {n_kept} cluster representatives")

        query_order, query_lengths = _read_query_ids_and_lengths(work_query)
        hits = run_blastn(work_query, db_prefix, threads=args.threads)
        classifications = classify_clusters(
            clusters, hits, query_lengths,
            n_subjects=len(subjects), total_analyzed_reads=total_analyzed_reads,
        )
        write_blast_report(clusters, classifications, args.out_dir, total_analyzed_reads)
    else:
        log("No --blast-subject provided; BLAST skipped (blast_result='.')")

    # Write top clusters with metadata + blast_result/flag columns.
    _write_top_clusters_tsv(
        clusters, args.out_dir, args.top,
        total_analyzed_reads,
        edit_dist=12 if loose_used else args.edit_dist,
        algorithm=args.algorithm,
        min_len=args.min_len,
        loose_criteria_used=loose_used,
        classifications=classifications,
    )
    log("done")


if __name__ == "__main__":
    main()
