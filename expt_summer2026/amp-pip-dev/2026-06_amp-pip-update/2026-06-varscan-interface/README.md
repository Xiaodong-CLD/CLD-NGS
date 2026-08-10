# 2026-06 — VarScan interface: `pileup2*` vs `mpileup2*`

**Type:** frozen investigation / decision record. Re-runs VarScan only; does not modify the pipeline.

## Purpose
Test whether `varscan_wrapper.py` can switch from the legacy `pileup2snp`/`pileup2indel` callers to `mpileup2snp`/`mpileup2indel` without changing calls.

**Decision:** Do NOT switch to `mpileup2*`. Indel calls differ in every sample/depth (missed positions, position shifts, VarFreq deltas up to ~4 pp). SNP calls differ at `d=8000`. The `mpileup2*` output format is also incompatible with all downstream CONS-format parsers in this repo.

## Required inputs
Two completed amp-pip-dev run directories:
- `<amp-pip-dev>/runs/T3_rebuild` — original run used for comparison.
- `<amp-pip-dev>/runs/T4_varscan_update` — second run; must contain per-sample `depth-sweep/` with pre-generated `d<D>.mpileup` files (produced by `../2026-06-mpileup-depth-cap/mpileup_depth_sweep.py`).

## Setup
```bash
conda activate variant-calling
```

## Run
```bash
cd expt_summer2026/amp-pip-dev/2026-06_amp-pip-update/2026-06-varscan-interface

# Fast (d=8000 only):
python compare_varscan_interface.py <amp-pip-dev>/runs/T3_rebuild \
                                    <amp-pip-dev>/runs/T4_varscan_update \
                                    --depths 8000

# Full (adds large ~1 GB d=0 mpileups; slow):
python compare_varscan_interface.py <amp-pip-dev>/runs/T3_rebuild \
                                    <amp-pip-dev>/runs/T4_varscan_update \
                                    --depths 0 8000
```

## Expected outputs
Per-sample/per-depth console report of interface-unique calls and VarFreq deltas for SNPs and indels separately. Raw per-call outputs cached under `runs/T4_varscan_update/`.

## Files
- `compare_varscan_interface.py` — interface-equivalence comparison.
