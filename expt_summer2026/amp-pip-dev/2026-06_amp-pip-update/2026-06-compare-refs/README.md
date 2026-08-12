# 2026-06 — insert-reference vs plasmid-reference variant concordance

**Type:** frozen investigation. Reads existing VarScan outputs only; does not modify the pipeline.

## Purpose
Sanity-check that insert-reference and plasmid-reference VarScan calls agree.
The insert reference is an exact sub-region of the plasmid (`plasmid_pos = insert_pos + 431`,
1-based, + strand). Verifies: no spurious backbone calls; concordance of calls inside the
insert window; VarFreq stability across the two references.

## Required inputs
Two completed pipeline run directories that share the same samples:
- `PLASMID_RESULTS_DIR` — e.g. `<amp-pip-dev>/runs/T1_5_threads/Results`
- `INSERT_RESULTS_DIR`  — e.g. `<amp-pip-dev>/runs/T_insertref/Results`

Each must contain `<sample>/varscan-output/<sample>.snp.vcf` and `*.indel.vcf`.

## Setup
```bash
conda activate variant-calling
```

## Run
```bash
cd expt_summer2026/amp-pip-dev/2026-06_amp-pip-update/2026-06-compare-refs
python compare_refs.py <amp-pip-dev>/runs/T1_5_threads/Results \
                       <amp-pip-dev>/runs/T_insertref/Results
```

## Expected outputs
Per-sample console report classifying each call as `MATCHED`, `PLASMID-ONLY`,
`INSERT-ONLY`, or `BACKBONE`; VarFreq drift flagged when Δ > 0.5 pp.

## Files
- `compare_refs.py` — insert↔plasmid call concordance report.
