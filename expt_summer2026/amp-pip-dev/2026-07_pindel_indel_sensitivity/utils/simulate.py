#!/usr/bin/env python3
"""Reference construction + read simulation for one benchmark condition.

This module owns the whole "produce the reads" front of the benchmark:
  1. carve a per-dataset amplicon window out of the WT reference,
  2. build a mutant amplicon carrying one known indel at a chosen breakpoint,
  3. ART-simulate paired-end amplicon reads from the WT and mutant references
     at fold counts matching a target allele fraction, and mix them into one
     sample (WT-only for AF=0 controls),
  4. write the run table (per-condition parameters + truth) consumed by
     utils/evaluate.py.

Reads are simulated at the post-trim length so fastp is skipped downstream; the
FASTQ naming (.R1/.R2.fastq.gz, no .trimmed) triggers the pipeline's no-trim
path. Coordinate convention: make_mutant_sequence returns the VCF-style anchor
(unchanged base immediately 5' of the event) so Pindel breakpoints line up.
"""

import csv
import gzip
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config  # noqa: E402


# Columns of the run table (results/manifest.csv), the single source of truth
# evaluate.py reads. Controls leave the indel columns blank.
RUN_TABLE_FIELDS = [
    "condition_id", "experiment", "dataset", "is_control", "indel_type",
    "length", "allele_freq", "pos", "read_len", "depth", "seed",
    "amplicon_pos", "rel_pos", "anchor_pos", "expected_svtype", "expected_svlen",
    "mutant_ref",
]


def run(cmd, **kwargs):
    """Run a shell command, raising on failure."""
    print(f"[run] {cmd}")
    subprocess.run(cmd, shell=True, check=True, **kwargs)


# --- FASTA helpers -----------------------------------------------------------

def read_single_fasta(path):
    """Return the concatenated, uppercased sequence of a single-record FASTA."""
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq).upper()


def write_fasta(path, name, sequence, width=70):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(sequence), width):
            fh.write(sequence[i:i + width] + "\n")


def extract_amplicon(dataset):
    """Return the WT amplicon sequence for a dataset.

    Reads the FASTA at ds["fasta"] and carves ds["region"] = (offset, length)
    when set, otherwise returns the full sequence.
    """
    ds = config.DATASETS[dataset]
    seq = read_single_fasta(ds["fasta"])
    region = ds.get("region")
    if region:
        offset, length = region
        if offset + length > len(seq):
            raise ValueError(
                f"Region [{offset}:{offset + length}] exceeds FASTA length "
                f"{len(seq)} for dataset '{dataset}'"
            )
        return seq[offset:offset + length]
    return seq


# --- mutant construction -----------------------------------------------------

def make_mutant_sequence(wt_amplicon, indel_type, length, pos):
    """Return (mutant_sequence, ref_allele, alt_allele, anchor_pos_1based).

    pos is the 1-based amplicon coordinate of the first affected base
    (deletion) or the insertion point. VCF-style anchor coordinates are
    returned for the truth table: the anchor is the unchanged base immediately
    5' of the event, so reported breakpoints from Pindel line up.
    """
    idx = pos - 1  # 0-based
    anchor_idx = idx - 1  # unchanged base before the event
    anchor_base = wt_amplicon[anchor_idx]
    anchor_pos_1based = anchor_idx + 1

    if indel_type == "del":
        deleted = wt_amplicon[idx:idx + length]
        mutant = wt_amplicon[:idx] + wt_amplicon[idx + length:]
        ref_allele = anchor_base + deleted
        alt_allele = anchor_base
    elif indel_type == "ins":
        inserted = config.INSERTION_TEMPLATE[:length]
        mutant = wt_amplicon[:idx] + inserted + wt_amplicon[idx:]
        ref_allele = anchor_base
        alt_allele = anchor_base + inserted
    else:
        raise ValueError(f"Unknown indel type: {indel_type}")

    return mutant, ref_allele, alt_allele, anchor_pos_1based


def _has_homopolymer(seq, center_idx, run_len=4, flank=3):
    """True if a homopolymer run of length >= run_len sits near center_idx."""
    lo = max(0, center_idx - flank)
    hi = min(len(seq), center_idx + flank + 1)
    window = seq[lo:hi]
    pat = r"A{%d,}|C{%d,}|G{%d,}|T{%d,}" % (run_len, run_len, run_len, run_len)
    return re.search(pat, window) is not None


def validate_positions(wt_amplicon, positions, label):
    """Warn if any 1-based breakpoint sits in/near a homopolymer run.

    Under exact-match scoring a breakpoint inside a repeat can be reported at
    the other edge of the homology and wrongly scored as a miss, so these are
    flagged for review. Returns the flagged positions.
    """
    flagged = [p for p in positions if _has_homopolymer(wt_amplicon, p - 1)]
    if flagged:
        print(f"[simulate] WARNING {label}: breakpoints in/near homopolymer "
              f"runs (ambiguous under exact scoring): {flagged}")
    return flagged


def dataset_amplicon_len(dataset):
    """Return the amplicon length for a dataset."""
    return len(extract_amplicon(dataset))


def ensure_wt_ref(dataset, outdir=config.OUTDIR):
    """Make sure a dataset's WT amplicon reference exists.

    Returns (seq, path, refs_dir).
    """
    refs_dir = os.path.join(outdir, "references")
    os.makedirs(refs_dir, exist_ok=True)
    wt_amplicon = extract_amplicon(dataset)
    stem = f"{config.amplicon_name(dataset)}_WT"
    wt_ref_path = os.path.join(refs_dir, f"{stem}.fa")
    if not os.path.exists(wt_ref_path):
        write_fasta(wt_ref_path, stem, wt_amplicon)
    return wt_amplicon, wt_ref_path, refs_dir


def build_mutant(wt_amplicon, dataset, indel_type, length, pos, refs_dir, outdir):
    """Build one mutant reference at a given position; return (path, truth_row).

    Reference filename includes dataset + position so distinct conditions get
    distinct references. Reuses an existing file when present.
    """
    mutant, ref_allele, alt_allele, anchor_pos = make_mutant_sequence(
        wt_amplicon, indel_type, length, pos
    )
    mut_name = f"{config.amplicon_name(dataset)}_{indel_type}_L{length}_p{pos}"
    mut_path = os.path.join(refs_dir, f"{mut_name}.fa")
    if not os.path.exists(mut_path):
        write_fasta(mut_path, mut_name, mutant)
    amplicon_len = len(wt_amplicon)
    truth_row = {
        "dataset": dataset,
        "indel_type": indel_type,
        "length": length,
        "expected_svtype": "DEL" if indel_type == "del" else "INS",
        "expected_svlen": -length if indel_type == "del" else length,
        "amplicon_pos": pos,            # first affected base (1-based)
        "amplicon_len": amplicon_len,   # for rel_pos computation in condition_row
        "anchor_pos": anchor_pos,       # VCF-style anchor (1-based)
        "ref_allele": ref_allele,
        "alt_allele": alt_allele,
        "mutant_ref": os.path.relpath(mut_path, outdir),
    }
    return mut_path, truth_row


# --- read simulation + mixing ------------------------------------------------

def art_simulate(ref_path, out_prefix, pairs, seed, read_id, read_len,
                 fragment_mean, fragment_sd):
    """Simulate standard paired-end reads from one reference.

    Uses ART -m/-s for randomly distributed fragments across the reference.
    `pairs` is the number of read pairs to generate.
    """
    if read_len >= fragment_mean:
        raise ValueError(
            f"ART requires read_len < fragment_mean, but "
            f"read_len={read_len} >= fragment_mean={fragment_mean}. "
            f"Exclude this read length from active_read_lengths."
        )
    cmd = (
        f"art_illumina -ss {config.SEQ_SYS} -p -na "
        f"-i {ref_path} -l {read_len} -c {pairs} "
        f"-m {fragment_mean} -s {fragment_sd} "
        f"-rs {seed} -d {read_id} -o {out_prefix} -q"
    )
    run(cmd)
    r1, r2 = f"{out_prefix}1.fq", f"{out_prefix}2.fq"
    for fq in (r1, r2):
        if not os.path.exists(fq):
            raise RuntimeError(
                f"ART did not produce expected output: {fq}\n"
                f"Check disk space (df -h) and the ART command above."
            )
    return r1, r2


def _concat_gzip(inputs, out_gz):
    with gzip.open(out_gz, "wt") as out:
        for path in inputs:
            with open(path) as fh:
                shutil.copyfileobj(fh, out)


def mix_counts(af, depth):
    """Split total read pairs into (wt_pairs, mut_pairs) for an allele fraction.

    mut_pairs = round(af * depth); wt_pairs is the remainder, so the two always
    sum to depth. Raises if the fraction rounds to zero mutant pairs.
    """
    mut_pairs = round(af * depth)
    if mut_pairs <= 0:
        raise ValueError(
            f"allele fraction {af} at depth {depth} yields 0 mutant pairs; "
            f"raise depth or AF."
        )
    return depth - mut_pairs, mut_pairs


def simulate_condition(cond_dir, mutant_ref, wt_ref, af, seed, depth, read_len,
                       fragment_mean, fragment_sd):
    """Produce mixed <cond>.R1.fastq.gz / .R2.fastq.gz for one condition."""
    cond_id = os.path.basename(cond_dir)
    wt_pairs, mut_pairs = mix_counts(af, depth)

    tmp = os.path.join(cond_dir, "art_tmp")
    os.makedirs(tmp, exist_ok=True)
    wt1, wt2 = art_simulate(wt_ref, os.path.join(tmp, "wt"), wt_pairs, seed,
                            "wt", read_len, fragment_mean, fragment_sd)
    mut1, mut2 = art_simulate(mutant_ref, os.path.join(tmp, "mut"), mut_pairs,
                              seed + 1, "mut", read_len, fragment_mean, fragment_sd)

    r1 = os.path.join(cond_dir, f"{cond_id}.R1.fastq.gz")
    r2 = os.path.join(cond_dir, f"{cond_id}.R2.fastq.gz")
    _concat_gzip([wt1, mut1], r1)
    _concat_gzip([wt2, mut2], r2)
    shutil.rmtree(tmp, ignore_errors=True)
    return r1, r2


def simulate_control(cond_dir, wt_ref, seed, depth, read_len,
                     fragment_mean, fragment_sd):
    """Produce WT-only <cond>.R1/.R2.fastq.gz for an AF=0 control condition."""
    cond_id = os.path.basename(cond_dir)
    tmp = os.path.join(cond_dir, "art_tmp")
    os.makedirs(tmp, exist_ok=True)
    wt1, wt2 = art_simulate(wt_ref, os.path.join(tmp, "wt"), depth, seed,
                            "wt", read_len, fragment_mean, fragment_sd)
    r1 = os.path.join(cond_dir, f"{cond_id}.R1.fastq.gz")
    r2 = os.path.join(cond_dir, f"{cond_id}.R2.fastq.gz")
    _concat_gzip([wt1], r1)
    _concat_gzip([wt2], r2)
    shutil.rmtree(tmp, ignore_errors=True)
    return r1, r2


# --- run table (manifest) ----------------------------------------------------

def condition_row(cond_id, experiment, truth_row, af, read_len, depth, seed):
    """Assemble a mutant run-table row from a build_mutant truth_row + params."""
    amplicon_len = truth_row.get("amplicon_len", 0)
    amplicon_pos = truth_row["amplicon_pos"]
    rel_pos = round(amplicon_pos / amplicon_len, 4) if amplicon_len else ""
    return {
        "condition_id": cond_id,
        "experiment": experiment,
        "dataset": truth_row["dataset"],
        "is_control": 0,
        "indel_type": truth_row["indel_type"],
        "length": truth_row["length"],
        "allele_freq": af,
        "pos": amplicon_pos,
        "read_len": read_len,
        "depth": depth,
        "seed": seed,
        "amplicon_pos": amplicon_pos,
        "rel_pos": rel_pos,
        "anchor_pos": truth_row["anchor_pos"],
        "expected_svtype": truth_row["expected_svtype"],
        "expected_svlen": truth_row["expected_svlen"],
        "mutant_ref": truth_row["mutant_ref"],
    }


def control_row(cond_id, experiment, dataset, read_len, depth, seed):
    """Assemble a WT-only (AF=0) control run-table row."""
    return {
        "condition_id": cond_id,
        "experiment": experiment,
        "dataset": dataset,
        "is_control": 1,
        "indel_type": "",
        "length": "",
        "allele_freq": 0.0,
        "pos": "",
        "read_len": read_len,
        "depth": depth,
        "seed": seed,
        "amplicon_pos": "",
        "rel_pos": "",
        "anchor_pos": "",
        "expected_svtype": "",
        "expected_svlen": "",
        "mutant_ref": "",
    }


def write_run_table(rows, outdir, experiment):
    """Write one experiment's run table (results/manifest.csv).

    Returns the path.
    """
    results_dir = os.path.join(outdir, "results")
    os.makedirs(results_dir, exist_ok=True)
    manifest_path = os.path.join(results_dir, "manifest.csv")
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=RUN_TABLE_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return manifest_path
