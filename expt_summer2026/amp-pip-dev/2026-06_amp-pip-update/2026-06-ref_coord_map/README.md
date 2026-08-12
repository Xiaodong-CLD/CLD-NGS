# 2026-06 — insert→plasmid coordinate map + CDS-annotated coverage

**Type:** frozen, dataset-specific investigation. Reads existing run outputs only; does not modify the pipeline.

## Purpose
Two tools for the 2026-06-09 P3 TLA dataset (`runs/P3_TLA`, GS/HC/LC plasmid):

1. `ref_coord_map.py` — maps 1-based insert-reference positions onto the full plasmid. Handles circular plasmids and falls back to end-anchor matching when the insert differs from the plasmid by edits/SNPs.
2. `analysis.py` — replots per-sample coverage from `<sample>.depth.txt` with HC/LC/GS CDS regions shaded, and writes `<sample>.cds_depth.txt`. Also used by `../2026-06-mpileup-depth-cap/mpileup_depth_sweep.py` (imports `read_consensus` and `pad_to_reference`).

## Required inputs
- FASTA files `PLASMID.fa` and `INSERT.fa` (for `ref_coord_map.py`).
- A completed amp-pip-dev run directory with `Results/<sample>/<sample>.depth.txt` and `reference/*.gbk` (for `analysis.py`). Original run: `<amp-pip-dev>/runs/P3_TLA`.

## Setup
```bash
conda activate variant-calling   # must include Biopython and plotly for analysis.py
```

## Run
```bash
cd expt_summer2026/amp-pip-dev/2026-06_amp-pip-update/2026-06-ref_coord_map

# Coordinate map (prints offset and strand to stdout):
python ref_coord_map.py <amp-pip-dev>/runs/P3_TLA/reference/PLASMID.fa \
                        <amp-pip-dev>/runs/P3_TLA/reference/INSERT.fa

# CDS-annotated coverage replot:
python analysis.py <amp-pip-dev>/runs/P3_TLA
```

## Expected outputs
- `ref_coord_map.py`: prints insert→plasmid offset and strand (or end-anchor positions) to stdout.
- `analysis.py`: writes `Results/<sample>/<sample>.cds_depth.txt` and `<sample>.cds_coverage.jpeg` inside the run directory.

## Files
- `ref_coord_map.py` — insert→plasmid coordinate mapper (circular-aware).
- `analysis.py` — CDS-shaded coverage replot + per-feature depth blocks; also imported by `../2026-06-mpileup-depth-cap/mpileup_depth_sweep.py`.
