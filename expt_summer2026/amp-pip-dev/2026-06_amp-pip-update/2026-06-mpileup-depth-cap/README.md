# 2026-06 — samtools mpileup `-d` depth cap vs VarScan output

**Type:** frozen investigation / decision record. Reads existing run outputs only; does not modify the pipeline.

## Purpose
Determine whether the samtools default `-d 8000` cap silently discards reads and biases VarScan variant calls on deep amplicon data. Tested values: `d=0` (unlimited), `8000`, `12000`, `16000` against run `T3_rebuild` (Lib1/Lib2/Lib3).

**Decision:** Use `-d 0` (uncapped) for this assay. The default `-d 8000` discards >60–77% of positions, inflates borderline ≥2% variant calls by 1–3 per sample, and never recovers from this bias at any finite cap tested.

## Required inputs
A completed amp-pip-dev pipeline run with sorted BAMs and a `reference/` directory:
- `<amp-pip-dev>/runs/T3_rebuild` — used for the original sweep.

## Setup
```bash
conda activate variant-calling
```
`mpileup_depth_sweep.py` also requires `plotly` (included in the `variant-calling` env).

## Run
```bash
cd expt_summer2026/amp-pip-dev/2026-06_amp-pip-update/2026-06-mpileup-depth-cap

# Step 1: generate sweep outputs (per-sample mpileup + pileup2cns for each -d value).
# Writes Results/<sample>/depth-sweep/ inside the run directory. Slow (re-runs pileup2cns).
python mpileup_depth_sweep.py <amp-pip-dev>/runs/T3_rebuild

# Step 2: compare depth ratios and VarFreq vs d=0 (stdlib only, fast).
python compare_d_modes.py <amp-pip-dev>/runs/T3_rebuild
```

## Expected outputs
- `Results/<sample>/depth-sweep/<sample>.d<D>.{mpileup,consensus.vcf,consensus_coverage.jpeg}` for each `-d` value — written inside the run directory.
- Console table: max consensus depth, fraction of positions downsampled, VarFreq deltas, and gained/lost calls per sample per `-d`.

## Files
- `mpileup_depth_sweep.py` — per-sample/per-`-d` mpileup + pileup2cns sweep. Imports `read_consensus` and `pad_to_reference` from `../2026-06-ref_coord_map/analysis.py`.
- `compare_d_modes.py` — depth-ratio / VarFreq / threshold comparison vs `d=0` (stdlib only).
- `compare_sweep_calls.py` — extended per-call variant comparison (requires pandas); optional.


