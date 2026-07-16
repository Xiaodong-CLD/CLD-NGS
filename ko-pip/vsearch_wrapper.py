#!/usr/bin/env python
"""VSEARCH consensus wrapper: paired reads -> top haplotypes -> whole-amplicon consensus.

Parametrized Python port of the legacy run_vsearch.sh (now in archive/),
followed by the haplotype caller in call_top_sequence_variants.py. For one
sample:

  1. vsearch --fastq_mergepairs   : merge R1/R2 into full-length amplicon reads.
  2. vsearch --fastx_uniques      : dereplicate exact merged reads (--sizeout).
  3. vsearch --sortbysize         : sort by abundance, drop singletons < min_size.
  4. top-N extraction             : keep the top_n most-abundant unique sequences.
  5. call_top_sequence_variants   : align each top sequence to the amplicon
                                    reference, group into whole-amplicon
                                    haplotypes (edge-masked), and write the
                                    consensus FASTA (multiple haplotypes, each
                                    with read support size + frequency).

Outputs:
  out_dir/<prefix>.consensus.fa                    (haplotype deliverable)
  scratch_dir/<prefix>.merged.fq,
  scratch_dir/<prefix>.derep.fa,
  scratch_dir/<prefix>.derep.sorted.fa,
  scratch_dir/<prefix>.top{N}.fa                   (vsearch intermediates)
  scratch_dir/<prefix>.variants.tsv, .summary.tsv, .alignments.txt

Requires: vsearch (conda), and call_top_sequence_variants.py alongside this file.

No-clobber: if <prefix>.consensus.fa already exists and is non-empty, the whole
stage is skipped (matches the other wrappers).
"""

import argparse
import os
import re
import subprocess

from pipeline_io import is_nonempty
from call_top_sequence_variants import main as call_top_sequence_variants


def _run(cmd):
    """Run a command list, echoing it with the pipeline prefix."""
    print("[Takeda Pipeline Message] " + " ".join(cmd))
    subprocess.check_call(cmd)


def _extract_top_n(sorted_fa, top_fa, top_n):
    """Copy the first ``top_n`` FASTA records (already size-sorted) into top_fa."""
    written = 0
    with open(sorted_fa) as fin, open(top_fa, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                if written >= top_n:
                    break
                written += 1
            fout.write(line)
    return written


def _sum_derep_sizes(derep_fa):
    """Sum the ``;size=N;`` counts across a dereplicated FASTA = total merged reads."""
    total = 0
    with open(derep_fa) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            m = re.search(r"size=(\d+)", line)
            if m:
                total += int(m.group(1))
    return total

def vsearch_runner(
    r1_path,
    r2_path,
    reference_path,
    out_dir,
    scratch_dir,
    sample_id,
    threads=8,
    top_n=20,
    min_size=20,
    consensus_min_size=50,
    consensus_min_freq=0.01,
    consensus_edge_bp=5,
):
    """Merge/derep/rank paired reads with VSEARCH, then call consensus haplotypes.

    Only the consensus FASTA deliverable lands in ``out_dir``; the vsearch
    intermediates and the per-query diagnostic tables are written to
    ``scratch_dir`` (tmp) with a ``<sample_id>.`` prefix. Returns the consensus
    FASTA path, or False on failure.
    """
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(scratch_dir, exist_ok=True)

    merged_fq = os.path.join(scratch_dir, f"{sample_id}.merged.fq")
    derep_fa = os.path.join(scratch_dir, f"{sample_id}.derep.fa")
    sorted_fa = os.path.join(scratch_dir, f"{sample_id}.derep.sorted.fa")
    top_fa = os.path.join(scratch_dir, f"{sample_id}.top{top_n}.fa")
    out_prefix = os.path.join(out_dir, sample_id)
    diag_prefix = os.path.join(scratch_dir, sample_id)
    consensus_fa = out_prefix + ".consensus.fa"

    # No-clobber: a complete consensus FASTA means the stage already ran.
    if is_nonempty(consensus_fa):
        print(f"[Takeda Pipeline Message] VSEARCH consensus already exists for {sample_id}; skipping.")
        return consensus_fa

    try:
        # Step 1: merge paired reads into full-length amplicon reads.
        print("[Takeda Pipeline Message] VSEARCH: merging paired reads...")
        _run([
            "vsearch",
            "--fastq_mergepairs", r1_path,
            "--reverse", r2_path,
            "--fastqout", merged_fq,
            "--fastq_minovlen", "20",
            "--fastq_maxdiffs", "5",
            "--fastq_maxee", "1.0",
            "--fastq_maxns", "0",
            "--threads", str(threads),
        ])

        # Step 2: dereplicate exact full-length merged reads.
        print("[Takeda Pipeline Message] VSEARCH: dereplicating merged reads...")
        _run([
            "vsearch",
            "--fastx_uniques", merged_fq,
            "--fastaout", derep_fa,
            "--sizeout",
            "--fasta_width", "0",
            "--threads", str(threads),
        ])

        # Step 3: sort by abundance, dropping rare uniques below min_size.
        print("[Takeda Pipeline Message] VSEARCH: sorting by abundance...")
        _run([
            "vsearch",
            "--sortbysize", derep_fa,
            "--output", sorted_fa,
            "--sizein",
            "--fasta_width", "0",
            "--minsize", str(min_size),
        ])

        # Step 4: keep the top-N most abundant sequences.
        n_top = _extract_top_n(sorted_fa, top_fa, top_n)
        print(f"[Takeda Pipeline Message] VSEARCH: extracted top {n_top} sequence(s) -> {top_fa}")
        if n_top == 0:
            print("[Takeda Pipeline Message] VSEARCH: no sequences passed min_size; no consensus produced.")
            return False

        # Total merged reads (frequency denominator) = summed derep sizes.
        total_reads = _sum_derep_sizes(derep_fa)

        # Step 5: build whole-amplicon consensus haplotypes from the top sequences.
        print("[Takeda Pipeline Message] Building consensus haplotypes from top sequences...")
        call_top_sequence_variants([
            "--ref", reference_path,
            "--query", top_fa,
            "--out_prefix", out_prefix,
            "--diag_prefix", diag_prefix,
            "--consensus_min_size", str(consensus_min_size),
            "--consensus_min_freq", str(consensus_min_freq),
            "--consensus_edge_bp", str(consensus_edge_bp),
            "--total_reads", str(total_reads),
        ])

        return consensus_fa

    except subprocess.CalledProcessError:
        print("[Takeda Pipeline Message] VSEARCH consensus processing failed.")
        return False


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Merge/derep/rank paired reads with VSEARCH and call whole-amplicon consensus haplotypes."
    )
    parser.add_argument("-1", "--r1", required=True, help="Forward (R1) FASTQ[.gz].")
    parser.add_argument("-2", "--r2", required=True, help="Reverse (R2) FASTQ[.gz].")
    parser.add_argument("-r", "--reference", required=True, help="Amplicon reference FASTA (single sequence).")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for the consensus FASTA deliverable.")
    parser.add_argument("--scratch_dir", required=True, help="Scratch directory (tmp) for vsearch intermediates and diagnostic tables.")
    parser.add_argument("-s", "--sample_id", required=True, help="Sample id (output file prefix).")
    parser.add_argument("--top_n", type=int, default=20, help="Number of top-abundance sequences to keep. Default 20.")
    parser.add_argument("--consensus_min_size", type=int, default=50, help="Minimum summed read support for a reported haplotype. Default 50.")
    parser.add_argument("--consensus_min_freq", type=float, default=0.01, help="Minimum frequency for a reported haplotype. Default 0.01 (1%%); set 0 to report all.")
    args = parser.parse_args(argv)

    result = vsearch_runner(
        r1_path=args.r1,
        r2_path=args.r2,
        reference_path=args.reference,
        out_dir=args.out_dir,
        scratch_dir=args.scratch_dir,
        sample_id=args.sample_id,
        top_n=args.top_n,
        consensus_min_size=args.consensus_min_size,
        consensus_min_freq=args.consensus_min_freq,
    )
    if not result:
        raise SystemExit(1)
    print(f"[Takeda Pipeline Message] VSEARCH consensus written: {result}")


if __name__ == "__main__":
    main()
