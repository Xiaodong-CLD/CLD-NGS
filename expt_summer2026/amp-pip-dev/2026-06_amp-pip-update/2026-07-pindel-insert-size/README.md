# 2026-07 — Pindel insert-size config (`398`) sensitivity sweep

**Type:** parameter-tuning investigation / decision record.

## Purpose
`pindel_wrapper.py` hard-codes the Pindel config insert-size field to `398`. This experiment sweeps that value over `[200, 300, 398, 500, 600]` plus each sample's empirical mean insert size (from `samtools stats`) and checks whether the choice affects the indel calls that reach the `*.100.vcf` summary (SI + LI + filtered D).

**Decision:** Keep `398` — it is arbitrary but harmless for this assay. Variant positions, ALT alleles, and GT are identical at all insert values tested; only sub-0.2% AD jitter is observed, well below the 2% cutoff. The empirical insert size for `T3_rebuild` is ~219–244, well below 398, confirming 398 was not derived from this library, but the field is insensitive to this difference for small indels on a short amplicon.

**Caveat:** insensitivity is specific to ~10 kb amplicon / ultra-deep / 2×75 / small SI. Not validated for large deletions/insertions.

## Required inputs
A completed amp-pip-dev run with sorted BAMs. Original run: `<amp-pip-dev>/runs/T3_rebuild`.

Reference FASTA is auto-selected from `<run_dir>/reference/*.fa`.

## Setup
```bash
conda activate variant-calling
```

## Run
```bash
cd expt_summer2026/amp-pip-dev/2026-06_amp-pip-update/2026-07-pindel-insert-size
python pindel_insert_sweep.py <amp-pip-dev>/runs/T3_rebuild
```

Compute-heavy: Pindel runs once per (sample × insert value). Resume-safe: skips Pindel if `_SI` file already exists.

## Expected outputs
- Per-sample `samtools stats` read-length and insert-size summary printed to stdout.
- `Results/<sample>/pindel-insert-sweep/insert-<V>/` for each value V — Pindel + pindel2vcf outputs written inside the run directory.
- Per-sample table of `insert_value → (SI, LI, D_filtered, total)` and called positions printed to stdout.
- Captured run output: `sweep_T3_rebuild.log`.

## Files
- `pindel_insert_sweep.py` — empirical insert-size measurement + Pindel sweep. Imports `filter_deletions_from_pindel2vcf` from `<amp-pip-dev>/filter_deletions.py`.
- `sweep_T3_rebuild.log` — captured stdout from the original sweep run.
