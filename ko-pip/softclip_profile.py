"""Standalone per-position soft-clip profile for an amplicon BAM.

Independent of `pipeline.py` / `downstream_analysis.py`. Produces a
coverage-style PNG track (`softclip_profile.png`) plotting, per amplicon
position, the number of reads whose alignment is soft-clipped on the left
(orange) and on the right (green). This module does NOT emit indel calls.

How to read the plot: a sharp local peak in left- or right-soft-clip reads
flags a candidate breakpoint — alignments terminate there because the true
sequence diverges (large insertion, large/complex deletion, structural edit).
Left and right peaks that flank a region bracket the edited/missing stretch.
Judge a peak by its height relative to local coverage (`softclip_fraction` in
the returned rows), not raw counts; a low, flat, evenly spread background at
the amplicon ends is ordinary primer/adapter trimming, not a KO signal.

External coordinates are 1-based; internal pysam coordinates are 0-based
half-open.
"""

from __future__ import annotations

import argparse
import os
import statistics
import sys
from typing import List, Optional, Tuple
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# BAM CIGAR op codes (SAM spec).
_OP_M, _OP_I, _OP_D, _OP_N, _OP_S, _OP_H, _OP_P, _OP_EQ, _OP_X = range(9)
_OPS_CONSUME_BOTH = {_OP_M, _OP_EQ, _OP_X}
_OPS_CONSUME_REF = {_OP_D, _OP_N}


def _read_usable(read, min_mapq: int) -> bool:
    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
        return False
    if (read.mapping_quality or 0) < min_mapq:
        return False
    if not read.cigartuples:
        return False
    return True


def _first_nonhard(cig):
    for op, ln in cig:
        if op != _OP_H:
            return op, ln
    return None, 0


def _last_nonhard(cig):
    for op, ln in reversed(cig):
        if op != _OP_H:
            return op, ln
    return None, 0


def _resolve_chrom(bamf: pysam.AlignmentFile, chrom: Optional[str]) -> str:
    refs = list(bamf.references)
    if chrom is not None:
        if chrom not in refs:
            raise ValueError(f"chrom {chrom!r} not in BAM header (have: {refs})")
        return chrom
    if len(refs) != 1:
        raise ValueError(
            f"BAM has {len(refs)} contigs ({refs}); pass --chrom explicitly."
        )
    return refs[0]


def build_softclip_profile(
    bam_path: str,
    output_dir: str,
    *,
    chrom: Optional[str] = None,
    sample_id: Optional[str] = None,
    min_mapq: int = 20,
    min_softclip_len: int = 10,
) -> Tuple[List[dict], Optional[str]]:
    """Scan a coordinate-sorted BAM and write a soft-clip profile PNG.

    Returns (per_position_rows, png_path_or_None).
    """
    os.makedirs(output_dir, exist_ok=True)
    if sample_id is None:
        sample_id = os.path.basename(bam_path).replace('.sorted.bam', '').replace('.bam', '')

    bamf = pysam.AlignmentFile(bam_path, 'rb')
    try:
        target = _resolve_chrom(bamf, chrom)
        ref_len = bamf.get_reference_length(target)

        m_depth = [0] * ref_len
        sc_left = [0] * ref_len
        sc_right = [0] * ref_len
        left_lens: List[List[int]] = [[] for _ in range(ref_len)]
        right_lens: List[List[int]] = [[] for _ in range(ref_len)]

        for read in bamf.fetch(target):
            if not _read_usable(read, min_mapq):
                continue
            cig = read.cigartuples
            ref_pos = read.reference_start

            for op, ln in cig:
                if op in _OPS_CONSUME_BOTH:
                    end = ref_pos + ln
                    if end > ref_len:
                        end = ref_len
                    for p in range(ref_pos, end):
                        m_depth[p] += 1
                    ref_pos += ln
                elif op in _OPS_CONSUME_REF:
                    ref_pos += ln
                # I/S/H/P do not advance ref_pos.

            first_op, first_ln = _first_nonhard(cig)
            if first_op == _OP_S and first_ln >= min_softclip_len:
                anchor = read.reference_start
                if 0 <= anchor < ref_len:
                    sc_left[anchor] += 1
                    left_lens[anchor].append(first_ln)

            last_op, last_ln = _last_nonhard(cig)
            if last_op == _OP_S and last_ln >= min_softclip_len:
                anchor = read.reference_end - 1
                if 0 <= anchor < ref_len:
                    sc_right[anchor] += 1
                    right_lens[anchor].append(last_ln)
    finally:
        bamf.close()

    sc_total = [sc_left[i] + sc_right[i] for i in range(ref_len)]
    sc_fraction = [
        (sc_total[i] / m_depth[i]) if m_depth[i] > 0 else 0.0
        for i in range(ref_len)
    ]
    med_left = [statistics.median(left_lens[i]) if left_lens[i] else 0 for i in range(ref_len)]
    med_right = [statistics.median(right_lens[i]) if right_lens[i] else 0 for i in range(ref_len)]

    # Per-position rows held in memory (not written to disk); the plot below
    # and any caller-side checks consume these directly.
    rows = [
        {
            'sample_id': sample_id, 'amplicon_id': target, 'position': i + 1,
            'm_depth': m_depth[i],
            'left_softclip_reads': sc_left[i],
            'right_softclip_reads': sc_right[i],
            'total_softclip_reads': sc_total[i],
            'softclip_fraction': sc_fraction[i],
            'median_left_softclip_len': med_left[i],
            'median_right_softclip_len': med_right[i],
        }
        for i in range(ref_len)
    ]

    png_path = os.path.join(output_dir, 'softclip_profile.png')
    try:
        positions = list(range(1, ref_len + 1))
        fig, ax = plt.subplots(1, 1, figsize=(10, 4))
        ax.plot(positions, sc_left, color='tab:orange', linewidth=1, label='left soft-clip')
        ax.plot(positions, sc_right, color='tab:green', linewidth=1, label='right soft-clip')
        ax.set_ylabel('soft-clip reads')
        ax.set_xlabel('amplicon position (1-based)')
        ax.set_title(f'{sample_id} | {target} | soft-clip profile')
        ax.legend(loc='best', fontsize=8)
        fig.tight_layout()
        fig.savefig(png_path, dpi=120)
        plt.close(fig)
        print(f'[Takeda Pipeline Message] Soft-clip profile plot: {png_path}')
    except Exception as exc:
        print(f'[Takeda Pipeline Message] Plot skipped ({exc}).')
        png_path = None

    return rows, png_path


def _parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--bam', help='Coordinate-sorted BAM with index.')
    p.add_argument('--outdir', help='Output directory.')
    p.add_argument('--chrom', default=None,
                   help='Amplicon contig name. Auto-detected if BAM has one contig.')
    p.add_argument('--sample-id', default=None, help='Sample id for the plot.')
    p.add_argument('--min-mapq', type=int, default=20)
    p.add_argument('--min-softclip-len', type=int, default=10)
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = _parse_args(argv)
    if not args.bam or not args.outdir:
        print('error: --bam and --outdir are required.', file=sys.stderr)
        return 2
    build_softclip_profile(
        bam_path=args.bam,
        output_dir=args.outdir,
        chrom=args.chrom,
        sample_id=args.sample_id,
        min_mapq=args.min_mapq,
        min_softclip_len=args.min_softclip_len,
    )
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
