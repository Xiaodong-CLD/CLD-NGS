# 2026-06 amp-pip-dev update experiments

Five frozen investigations from June–July 2026 into the amp-pip-dev amplicon
pipeline. Each experiment is self-contained under its own subdirectory with a
README, scripts, and (where applicable) a captured log.

## Experiments

| Directory | Topic |
|---|---|
| [2026-06-compare-refs](2026-06-compare-refs/README.md) | Insert-ref vs plasmid-ref VarScan call concordance |
| [2026-06-mpileup-depth-cap](2026-06-mpileup-depth-cap/README.md) | samtools `mpileup -d` cap effect on VarScan calls |
| [2026-06-ref_coord_map](2026-06-ref_coord_map/README.md) | Insert→plasmid coordinate map + CDS-annotated coverage |
| [2026-06-varscan-interface](2026-06-varscan-interface/README.md) | `pileup2*` vs `mpileup2*` VarScan interface equivalence |
| [2026-07-pindel-insert-size](2026-07-pindel-insert-size/README.md) | Pindel insert-size config (`398`) sensitivity sweep |

Pipeline source code lives in `ngs/amp-pip-dev/` and is not copied here.

---

## Testing

### What the tests cover

`tests/test_experiment_structure.py` — static checks (no execution):
- All experiment directories, README files, and scripts exist.
- All Python scripts parse without syntax errors.
- Cross-experiment import paths resolve correctly:
  - `mpileup_depth_sweep.py` → `expt_summer2026/amp-pip-dev/2026-06_amp-pip-update/2026-06-ref_coord_map/analysis.py`
  - `pindel_insert_sweep.py` → `amp-pip-dev/filter_deletions.py`
- Each README contains a `## Run` and `## Required inputs` section.

`tests/test_experiment_smoke.py` — minimal functional tests using synthetic data:
- `compare_refs.py` — full end-to-end run with synthetic VarScan tables.
- `compare_d_modes.py` — full end-to-end run with synthetic consensus VCFs (stdlib only).
- `compare_sweep_calls.py` — full run writing output TSVs; column names validated (requires pandas).
- `ref_coord_map.py` — exact-match and anchor-fallback cases on tiny synthetic FASTAs.
- `compare_varscan_interface.py` — `--help` exit check (argparse setup).
- `mpileup_depth_sweep.py` — import-path check only (samtools + varscan needed for full run).
- `pindel_insert_sweep.py` — import-path check only (pindel + samtools needed for full run).

### What the tests do NOT validate

- Scientific correctness of any result.
- Full pipeline runs on real datasets.
- `analysis.py` end-to-end (requires Biopython + plotly + a real GenBank reference).
- `mpileup_depth_sweep.py` end-to-end (requires samtools, varscan, sorted BAMs).
- `pindel_insert_sweep.py` end-to-end (requires pindel, samtools, sorted BAMs).
- `compare_varscan_interface.py` end-to-end (requires varscan and pre-generated mpileup files).

Tests that need absent tools are skipped automatically with a message naming the missing tool.

### Requirements

`pytest` and `pandas` must be installed in the active environment. The tests are designed to run under the `variant-calling` conda environment:

```bash
conda activate variant-calling
conda install pytest pandas
```

### How to run

```bash
cd ngs/expt_summer2026/amp-pip-dev/2026-06_amp-pip-update

# Run the full suite:
bash run_tests.sh

# Retain temporary outputs when a test fails (useful for debugging):
bash run_tests.sh --keep-failed
```

### Optional full integration tests

The end-to-end tests for experiments that call samtools/varscan/pindel are not
part of the default suite.  To run them manually:

```bash
# mpileup depth sweep (needs variant-calling conda env + run T3_rebuild):
conda activate variant-calling
cd ngs/expt_summer2026/amp-pip-dev/2026-06_amp-pip-update/2026-06-mpileup-depth-cap
python mpileup_depth_sweep.py ngs/amp-pip-dev/runs/T3_rebuild
python compare_d_modes.py     ngs/amp-pip-dev/runs/T3_rebuild

# Pindel insert-size sweep (needs variant-calling env + run T3_rebuild):
conda activate variant-calling
cd ngs/expt_summer2026/amp-pip-dev/2026-06_amp-pip-update/2026-07-pindel-insert-size
python pindel_insert_sweep.py ngs/amp-pip-dev/runs/T3_rebuild

# VarScan interface comparison (needs variant-calling env + runs T3_rebuild + T4_varscan_update):
conda activate variant-calling
cd ngs/expt_summer2026/amp-pip-dev/2026-06_amp-pip-update/2026-06-varscan-interface
python compare_varscan_interface.py ngs/amp-pip-dev/runs/T3_rebuild \
                                    ngs/amp-pip-dev/runs/T4_varscan_update \
                                    --depths 8000
```
