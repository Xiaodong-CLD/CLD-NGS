# 2026-07 Previous Runs Comparison

## Goal

Verify that the upgraded amp-pip-dev pipeline produces equivalent variant calls
to the legacy MiSeq Analysis Server results for the same FASTQ inputs.

## Approach

Run the new pipeline locally on four datasets (AvB1_HC, AvB1_LC, TnT_CCR8_HC,
TnT_CCR8_LC), then compare VarScan SNP/indel calls position-by-position against
the old Results directories on the Z: drive. Report matched, gained, and lost
variants plus VarFreq deltas.

## How to reproduce

```bash
cd /home/robolab/cindy_ren/ngs/expt_summer2026/amp-pip-dev/2026-07_prev_runs_compare

# 1. Re-run the pipeline on each dataset to produce local Results/ dirs
#    (requires Z: drive mounted and variant-calling conda env active).

# 2. Compare old vs new
python compare_runs.py \
    --old-results "/mnt/z/CLD/MiSeq/MiSeq CLD NGS results/250409_M06815_0016_000000000-LMG4G/20250410_114445/Results-1" \
    --new-results AvB1_HC/Results/ --output-dir comparison/AvB1_HC/

python compare_runs.py \
    --old-results "/mnt/z/CLD/MiSeq/MiSeq CLD NGS results/250411_M06815_0017_000000000-L85D2/20250412_045938/Results" \
    --new-results AvB1_LC/Results/ --output-dir comparison/AvB1_LC/

python compare_runs.py \
    --old-results "/mnt/z/CLD/MiSeq/251117_M06815_0020_000000000-LWCMF/Alignment_1/Results-TnT_CCR8-HC" \
    --new-results TnT_CCR8_HC/Results/ --output-dir comparison/TnT_CCR8_HC/

python compare_runs.py \
    --old-results "/mnt/z/CLD/MiSeq/251117_M06815_0020_000000000-LWCMF/Alignment_1/Results-TnT_CCR8-LC" \
    --new-results TnT_CCR8_LC/Results/ --output-dir comparison/TnT_CCR8_LC/

# 3. Summarize
python comparison_summary.py comparison/AvB1_HC/
python comparison_summary.py comparison/AvB1_LC/
python comparison_summary.py comparison/TnT_CCR8_HC/
python comparison_summary.py comparison/TnT_CCR8_LC/

# 4. Run unit tests (no external data needed)
python test_comparison_scripts.py
```

## Scripts

| File | Purpose |
|------|---------|
| `compare_runs.py` | Per-sample SNP/indel comparison between two Results dirs |
| `comparison_summary.py` | Human-readable summary of a comparison output dir |
| `compare_depth_sweep.py` | Compare VarScan consensus at two mpileup -d caps |
| `winpath.py` | Convert Windows paths to WSL mount paths |
| `test_comparison_scripts.py` | Synthetic-data smoke tests for all scripts |
