#!/usr/bin/env python3
"""Extract, cluster, and BLAST soft-clipped sequences.

Stages: BAM walk -> starcode clustering (per 5'/3' side) -> BLAST top-N
representatives -> classified top_clusters.tsv.

Key output columns in top_clusters.tsv:
  blast_result     'hit' / 'partial' / 'no hit' (pident>=90%, cov>=80% for hit)
  construct_blast  hit / partial / '.' vs plasmid/construct subjects
  reference_blast  hit / partial / '.' vs amplicon reference (a ref match
                   suggests large deletion -- mapper couldn't place the clip)
  large_insert     'yes' if clip is long (>= MIN_FLAG_CLIP_LEN) and well
                   supported (>= MIN_SUPPORT_FRAC of total analyzed reads)

No-clobber: if top_clusters.tsv already exists the whole stage is skipped.
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

from pipeline_io import log, LOG_PREFIX, is_nonempty

EXCLUDE_FLAGS = 4 | 256 | 1024 | 2048  # 3332: unmapped/secondary/dup/supp
CIGAR_S = 4
BLAST_DB_EXTS = (".nhr", ".nin", ".nsq")

# Insert-flag thresholds (see classify_clusters).
FULL_HIT_MIN_PIDENT = 90.0   # percent identity for a "hit" (construct present)
FULL_HIT_MIN_COV = 80.0      # percent of the clip that must align for a "hit"
MIN_SUPPORT_FRAC = 0.01      # cluster reads / total analyzed reads (full BAM)
MIN_FLAG_CLIP_LEN = 30       # minimum representative clip length (bp) to flag
LARGE_DEL_ANCHOR_TOL_BP = 5  # bp: flag large_deletion if BLAST hit endpoint within this of anchor_pos

# Clustering / output hyperparameters.
# These are fixed algorithmic choices, not per-run tunable knobs.  Change them
# here if future validation motivates different defaults.
#
#   CLUSTER_TOP_N     = 20   Top-N clusters written to FASTA and BLAST-queried.
#                            20 captures the dominant junctions in practice while
#                            keeping the BLAST query set small.
#   CLUSTER_EDIT_DIST = 8    Levenshtein distance for starcode sphere-clustering.
#                            Chosen to be aggressive enough to collapse sequencing-
#                            error variants of the same junction into one cluster,
#                            without merging biologically distinct junctions
#                            (tested on Fut8 KO libraries; typical clip length
#                            30-170 bp, so d=8 is ~5-25% of clip length).
#   CLUSTER_ALGORITHM = sphere  Sphere-clustering is the most aggressive starcode
#                            mode: every read within edit_dist of a centroid joins
#                            that cluster, which gives the cleanest locus-level
#                            counts for amplicon soft-clip data.
#   CLUSTER_MIN_LEN   = 10   Minimum soft-clip length accepted into clustering.
#                            Clips shorter than this are overwhelmingly random
#                            mapping artefacts; 10 bp is the starcode minimum.
#   CLUSTER_THREADS   = 4    CPU threads for starcode and blastn.
CLUSTER_TOP_N = 20
CLUSTER_EDIT_DIST = 8
CLUSTER_ALGORITHM = "sphere"
CLUSTER_MIN_LEN = 10
CLUSTER_THREADS = 4


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
    """Run starcode; yield (centroid, total_count, members) per cluster."""
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
                (_, top_anchor), loc_n = cluster_anchors.most_common(1)[0]
                out.append({
                    "side": side,
                    "reads": total,
                    "members": len(members),
                    "rep_seq": centroid,
                    "rep_reads": sc.get(centroid, 0),
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
    """Write top_clusters.tsv with a comment-line metadata header and top-N rows.

    support_frac = cluster reads / total_reads (all BAM reads, post flag-filter).
    Classification columns (blast_result, construct_blast, reference_blast,
    large_insert) are filled when ``classifications`` is provided, else '.'.
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
        fh.write(f"# min_support_frac: {MIN_SUPPORT_FRAC:.6f}\n")
        # Write header (support_frac = cluster_reads / total_analyzed_reads)
        fh.write("rank\tcluster_reads\tsupport_frac\tside\tanchor\tlength\trep_sequence\tblast_result\tconstruct_blast\treference_blast\tflag\tvariant\n")
        # Write top-N clusters
        for rank, c in enumerate(clusters[:top_n], start=1):
            support_frac = c["reads"] / denom
            cls = classifications.get(rank, {})
            blast_result = cls.get("blast_result", ".")
            construct_blast = cls.get("construct_blast", ".")
            reference_blast = cls.get("reference_blast", ".")
            flag = cls.get("flag", ".")
            side = c["side"]
            anchor_int = int(c["top_anchor"])
            top_hits = cls.get("top_hits") or []
            if "large_deletion" in flag and top_hits:
                h = top_hits[0]
                vs = h["send"] + 1 if side == "5p" else anchor_int + 1
                ve = anchor_int - 1 if side == "5p" else h["sstart"] - 1
                variant = f"DEL;{vs},{ve};{ve - vs + 1}"
            elif "large_insert" in flag:
                variant = f"INS;{anchor_int},{anchor_int};{len(c['rep_seq'])}"
            else:
                variant = "."
            fh.write(
                f"{rank}\t{c['reads']}\t{support_frac:.6f}\t{side}\t"
                f"{anchor_int}\t{len(c['rep_seq'])}\t"
                f"{c['rep_seq']}\t{blast_result}\t"
                f"{construct_blast}\t{reference_blast}\t"
                f"{flag}\t{variant}\n"
            )
    log(f"Wrote {tsv_path} ({min(top_n, len(clusters))} clusters, {len(clusters)} total)")
    return tsv_path



def write_cluster_table_and_fasta(
    clusters: List[Dict],
    out_dir: pathlib.Path,
    top_n: int,
) -> pathlib.Path:
    """Write top_softclip_sequence.fa with the top-N cluster representatives (BLAST query)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = out_dir / "top_softclip_sequence.fa"
    total = sum(c["reads"] for c in clusters)
    denom = total if total else 1

    with fasta_path.open("w") as fh:
        for rank, c in enumerate(clusters[: top_n], start=1):
            freq = c["reads"] / denom
            # The sequence id (cluster{rank}) stays a bare first token so BLAST
            # keeps it as qseqid; the metadata fields are ';'-joined.
            header = (
                f">cluster{rank} side={c['side']};cluster_reads={c['reads']};"
                f"freq_frac={freq:.6f};len={len(c['rep_seq'])};"
                f"top_locus={c['top_anchor']};top_locus_reads={c['locus_reads']}"
            )
            fh.write(f"{header}\n{c['rep_seq']}\n")

    log(f"Wrote {fasta_path} (top {min(top_n, len(clusters))})")
    return fasta_path


def write_per_clip_fasta_and_metadata(
    bam_path: pathlib.Path,
    fasta_path: pathlib.Path,
    meta_path: pathlib.Path,
    min_len: int,
) -> int:
    """Write per-read soft-clip FASTA with unique query IDs and a metadata TSV.

    query_id = f"{read.query_name}|{side}|{index}"
    anchor_pos: 1-based ref pos of first aligned base (5p) or last aligned base (3p).
    Returns number of clips written.
    """
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    idx = 0
    with (fasta_path.open("w") as fa,
          meta_path.open("w") as mt,
          pysam.AlignmentFile(str(bam_path), "rb") as bam):
        mt.write("query_id\tanchor_chrom\tanchor_pos\tclip_side\n")
        for read in bam.fetch(until_eof=True):
            if read.flag & EXCLUDE_FLAGS:
                continue
            cigar = read.cigartuples
            seq = read.query_sequence
            if not cigar or not seq:
                continue
            chrom = read.reference_name or "*"
            # 5' left clip: anchor = first aligned reference base (1-based)
            if cigar[0][0] == CIGAR_S and cigar[0][1] >= min_len:
                clip_len = cigar[0][1]
                anchor = (read.reference_start + 1) if read.reference_start is not None else 0
                qid = f"{read.query_name}|5p|{idx}"
                fa.write(f">{qid}\n{seq[:clip_len]}\n")
                mt.write(f"{qid}\t{chrom}\t{anchor}\t5p\n")
                idx += 1
            # 3' right clip: anchor = last aligned reference base (1-based)
            if cigar[-1][0] == CIGAR_S and cigar[-1][1] >= min_len:
                clip_len = cigar[-1][1]
                anchor = read.reference_end if read.reference_end is not None else 0
                qid = f"{read.query_name}|3p|{idx}"
                fa.write(f">{qid}\n{seq[-clip_len:]}\n")
                mt.write(f"{qid}\t{chrom}\t{anchor}\t3p\n")
                idx += 1
    log(f"Wrote {idx} per-clip entries -> {fasta_path}, metadata -> {meta_path}")
    return idx


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
    ref_subject_ids: set | None = None,
) -> Dict[int, dict]:
    """Classify each BLASTed top cluster into blast_result / construct_blast /
    reference_blast / flag columns. Only clusters present in
    query_lengths (the BLASTed top-N) are classified; the rest are skipped.
    """
    from collections import defaultdict
    by_query: Dict[str, List[dict]] = defaultdict(list)
    for h in hits:
        by_query[h["qseqid"]].append(h)

    ref_ids = ref_subject_ids or set()
    denom = total_analyzed_reads if total_analyzed_reads else 1
    out: Dict[int, dict] = {}
    for idx, c in enumerate(clusters, start=1):
        qid = f"cluster{idx}"
        if qid not in query_lengths:
            continue  # only the BLASTed top-N clusters are classified/reported
        qlen = query_lengths.get(qid) or len(c["rep_seq"])

        best_per_subj: Dict[str, dict] = {}
        for h in by_query.get(qid, []):
            s = h["sseqid"]
            if s not in best_per_subj or h["bitscore"] > best_per_subj[s]["bitscore"]:
                best_per_subj[s] = h
        top_hits = sorted(best_per_subj.values(), key=lambda r: -r["bitscore"])[:max(1, n_subjects)]

        def _hit_level(hit_list):
            if not hit_list:
                return "no hit"
            is_full = any(
                h["pident"] >= FULL_HIT_MIN_PIDENT
                and (100.0 * h["length"] / qlen if qlen else 0.0) >= FULL_HIT_MIN_COV
                for h in hit_list
            )
            return "hit" if is_full else "partial"

        # Split hits into construct (plasmid) vs reference groups.
        construct_hits = [h for h in top_hits if h["sseqid"] not in ref_ids]
        ref_hits = [h for h in top_hits if h["sseqid"] in ref_ids]

        blast_result = _hit_level(top_hits)

        # construct_blast: does the clip match the known construct/plasmid?
        level = _hit_level(construct_hits)
        if level == "hit":
            construct_blast = "construct_match"
        elif level == "partial":
            construct_blast = "partial"
        else:
            construct_blast = "."

        # reference_blast: does the clip match the reference amplicon?
        # A match here suggests a large deletion (the soft-clipped bases are
        # reference-sequence that the aligner couldn't place in context).
        if ref_ids:
            ref_level = _hit_level(ref_hits)
            if ref_level == "hit":
                reference_blast = "reference_match"
            elif ref_level == "partial":
                reference_blast = "partial"
            else:
                reference_blast = "."
        else:
            reference_blast = "."

        # large_insert: long, well-supported, non-reference matching clip (independent of BLAST).
        support_frac = c["reads"] / denom
        large_insert = (
            support_frac >= MIN_SUPPORT_FRAC
            and len(c["rep_seq"]) >= MIN_FLAG_CLIP_LEN
            and reference_blast == "."
        )

        # large_deletion: well-supported clip whose top BLAST hit endpoint lands
        # LARGE_DEL_ANCHOR_TOL_BP+ away from the genomic anchor, and is a
        # reference match. (5p: compare send; 3p: compare sstart)
        anchor_pos_v = int(c.get("top_anchor") or 0)
        side_v = c.get("side", "5p")
        large_deletion = (
            support_frac >= MIN_SUPPORT_FRAC
            and reference_blast != "."
            and any(
                (abs(h["send"] - anchor_pos_v) if side_v == "5p" else abs(h["sstart"] - anchor_pos_v))
                >= LARGE_DEL_ANCHOR_TOL_BP
                for h in top_hits
            )
        )

        # flag: merge large_insert / large_deletion into a single column.
        # Both flags should not be reached under current logic.
        flag_parts = []
        if large_insert:
            flag_parts.append("large_insert")
        if large_deletion:
            flag_parts.append("large_deletion")
        flag = ", ".join(flag_parts) if flag_parts else "."

        out[idx] = {"blast_result": blast_result,
                    "construct_blast": construct_blast,
                    "reference_blast": reference_blast,
                    "flag": flag,
                    "top_hits": top_hits, "qlen": qlen}
    return out


def write_blast_report(
    clusters: List[Dict],
    classifications: Dict[int, dict],
    out_dir: pathlib.Path,
    total_reads: int,
) -> pathlib.Path:
    """Write blast_hits_full.tsv: one row per top-N hit per cluster (always written,
    including 'no hit' clusters)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    blast_tsv = out_dir / "blast_hits_full.tsv"
    header = [
        "cluster_id", "n_reads", "support_frac", "repseq_len",
        "blast_result", "blast_hit_sseqid", "pident", "evalue",
        "aln_len", "cov_pct", "sstart", "send",
        "anchor_pos", "clip_side", "large_deletion",
    ]
    denom = total_reads if total_reads else 1

    with blast_tsv.open("w") as fh:
        fh.write("\t".join(header) + "\n")
        for idx, c in enumerate(clusters, start=1):
            cls = classifications.get(idx)
            if cls is None:
                continue
            qid = f"cluster{idx}"
            support_frac = c["reads"] / denom
            qlen = cls["qlen"] or len(c["rep_seq"])
            repseq_len = len(c["rep_seq"])
            anchor_pos = c.get("top_anchor", ".")
            clip_side = c.get("side", ".")
            ap_int = int(anchor_pos) if anchor_pos != "." else 0
            if cls["top_hits"]:
                for h in cls["top_hits"]:
                    cov = 100.0 * h["length"] / qlen if qlen else 0.0
                    dist = abs(h["send"] - ap_int) if clip_side == "5p" else abs(h["sstart"] - ap_int)
                    row_ld = "yes" if dist <= LARGE_DEL_ANCHOR_TOL_BP else "no"
                    fh.write("\t".join([
                        qid, str(c["reads"]), f"{support_frac:.6f}", str(repseq_len),
                        cls["blast_result"], h["sseqid"], f"{h['pident']:.2f}",
                        f"{h['evalue']:.2e}", str(h["length"]), f"{cov:.1f}",
                        str(h["sstart"]), str(h["send"]),
                        str(anchor_pos), clip_side, row_ld,
                    ]) + "\n")
            else:
                fh.write("\t".join([
                    qid, str(c["reads"]), f"{support_frac:.6f}", str(repseq_len),
                    cls["blast_result"], ".", ".", ".", ".", ".",
                    ".", ".",
                    str(anchor_pos), clip_side, ".",
                ]) + "\n")
    log(f"Wrote {blast_tsv} ({len(classifications)} clusters classified)")
    return blast_tsv


def join_blast_to_metadata(
    hits: List[dict],
    meta_path: pathlib.Path,
    out_path: pathlib.Path,
) -> pathlib.Path:
    """Join per-clip BLAST hits to per_clip_metadata.tsv via qseqid; write TSV."""
    meta: Dict[str, Tuple[str, str, str]] = {}
    with meta_path.open() as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                meta[parts[0]] = (parts[1], parts[2], parts[3])
    out_cols = BLAST_FIELDS + ["anchor_chrom", "anchor_pos", "clip_side"]
    with out_path.open("w") as fh:
        fh.write("\t".join(out_cols) + "\n")
        for h in hits:
            qid = h["qseqid"]
            anchor_chrom, anchor_pos, clip_side = meta.get(qid, (".", ".", "."))
            row = [str(h.get(c, ".")) for c in BLAST_FIELDS]
            row += [anchor_chrom, anchor_pos, clip_side]
            fh.write("\t".join(row) + "\n")
    log(f"Wrote {len(hits)} per-clip BLAST rows -> {out_path}")
    return out_path


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        description="Extract, cluster, and BLAST soft-clipped sequences."
    )
    ap.add_argument("-b", "--bam", required=True, type=pathlib.Path,
                    help="Input BAM (uninformative or full sample).")
    ap.add_argument("-o", "--out-dir", required=True, type=pathlib.Path,
                    help="Output directory for kept deliverables (top_clusters.tsv, "
                         "top_softclip_sequence.fa).")
    ap.add_argument("--scratch-dir", type=pathlib.Path, default=None,
                    help="Directory for diagnostic BLAST outputs (blast_hits_full.tsv); "
                         "defaults to <out-dir>.")

    # BLAST knobs
    ap.add_argument("--blast-subject", default=None,
                    help="Local subject FASTA (comma-separated list) or omit to skip BLAST.")
    ap.add_argument("--ref-fasta", type=pathlib.Path, default=None,
                    help="Amplicon reference FASTA. When provided it is added to the BLAST "
                         "database so clips matching the reference populate reference_blast "
                         "(suggesting large deletions); construct subjects and reference "
                         "are classified separately.")
    ap.add_argument("--blast-db-dir", type=pathlib.Path, default=None,
                    help="Directory for BLAST DB cache (default: <out-dir>/blast_db).")
    args = ap.parse_args(argv)

    if not args.bam.exists():
        sys.exit(f"{LOG_PREFIX} ERROR: BAM not found: {args.bam}")

    # No-clobber: top_clusters.tsv is written last (after BLAST), so its
    # presence means the full stage already completed.
    top_clusters_path = args.out_dir / "top_clusters.tsv"
    if is_nonempty(top_clusters_path):
        log(f"top_clusters.tsv already exists; skipping starcode + BLAST.")
        return

    # ---- STAGE 1: Extract and cluster ----
    log(f"Scanning soft clips in {args.bam} "
        f"(min_len={CLUSTER_MIN_LEN}, edit_dist={CLUSTER_EDIT_DIST}, algo={CLUSTER_ALGORITHM})")
    records: List[Tuple[str, str, str, int]] = list(iter_softclips(args.bam, CLUSTER_MIN_LEN))
    total = len(records)
    log(f"Collected {total} soft-clip segments; clustering with starcode…")

    clusters = cluster_with_starcode(
        records, edit_dist=CLUSTER_EDIT_DIST, algorithm=CLUSTER_ALGORITHM, threads=CLUSTER_THREADS,
    )
    log(f"Formed {len(clusters)} clusters")
    
    # Retry with looser criteria if no clusters found
    loose_used = False
    if len(clusters) == 0 and len(records) > 0:
        log(f"No clusters with edit_dist={CLUSTER_EDIT_DIST}; retrying with edit_dist=12...")
        clusters = cluster_with_starcode(
            records, edit_dist=12, algorithm=CLUSTER_ALGORITHM, threads=CLUSTER_THREADS,
        )
        log(f"Formed {len(clusters)} clusters with loose criteria")
        loose_used = True
    
    if len(clusters) == 0:
        log("No clusters found even with loose criteria; skipping output files")
        return

    total_analyzed_reads = count_total_reads(args.bam)
    log(f"Total reads (post-filter) for support denominator = {total_analyzed_reads}")

    # Generate FASTA of cluster representatives (also the BLAST queries).
    fasta_path = write_cluster_table_and_fasta(clusters, args.out_dir, CLUSTER_TOP_N)
    scratch_dir = args.scratch_dir or args.out_dir

    # ---- STAGE 2: BLAST + insert classification ----
    subjects = _parse_subjects(args.blast_subject or "")
    ref_fasta = args.ref_fasta  # may be None
    classifications: Dict[int, dict] = {}
    if subjects or ref_fasta:
        all_subjects = list(subjects)
        ref_subject_ids: set = set()
        if ref_fasta:
            if not ref_fasta.exists():
                sys.exit(f"{LOG_PREFIX} ERROR: --ref-fasta not found: {ref_fasta}")
            all_subjects.append(ref_fasta)
            # Collect the sequence ID(s) from the reference FASTA so we can
            # separate reference hits from construct hits during classification.
            _, ref_lens = _read_query_ids_and_lengths(ref_fasta)
            ref_subject_ids = set(ref_lens.keys())
        blast_db_dir = args.blast_db_dir or (args.out_dir / "blast_db")
        db_prefix = prepare_blast_db(all_subjects, blast_db_dir)

        # BLAST the top-cluster representative FASTA directly.
        query_order, query_lengths = _read_query_ids_and_lengths(fasta_path)
        hits = run_blastn(fasta_path, db_prefix, threads=CLUSTER_THREADS)
        classifications = classify_clusters(
            clusters, hits, query_lengths,
            n_subjects=len(all_subjects), total_analyzed_reads=total_analyzed_reads,
            ref_subject_ids=ref_subject_ids,
        )
        write_blast_report(clusters, classifications, scratch_dir, total_analyzed_reads)
        # Per-read soft-clip FASTA + anchor-aware BLAST
        per_clip_fa = scratch_dir / "per_clip.fa"
        per_clip_meta = scratch_dir / "per_clip_metadata.tsv"
        n_clips = write_per_clip_fasta_and_metadata(
            args.bam, per_clip_fa, per_clip_meta, CLUSTER_MIN_LEN
        )
        if n_clips > 0:
            clip_hits = run_blastn(per_clip_fa, db_prefix, threads=CLUSTER_THREADS)
            join_blast_to_metadata(
                clip_hits, per_clip_meta, scratch_dir / "per_clip_blast.tsv"
            )
    else:
        log("No --blast-subject or --ref-fasta provided; BLAST skipped (blast_result='.')")

    # Write top clusters with metadata + classification columns.
    _write_top_clusters_tsv(
        clusters, args.out_dir, CLUSTER_TOP_N,
        total_analyzed_reads,
        edit_dist=12 if loose_used else CLUSTER_EDIT_DIST,
        algorithm=CLUSTER_ALGORITHM,
        min_len=CLUSTER_MIN_LEN,
        loose_criteria_used=loose_used,
        classifications=classifications,
    )
    log("done")


if __name__ == "__main__":
    main()
