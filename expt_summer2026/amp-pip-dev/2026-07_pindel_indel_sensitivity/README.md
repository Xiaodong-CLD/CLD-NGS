# Pindel indel-sensitivity benchmark

Simulation study measuring how well **Pindel** and other short-read indel
callers detect variants in amplicon-seq data.

## Key questions

| | Question | Script |
|---|---|---|
| **A** | At what allele fraction and indel length does Pindel start to miss variants? | `detect_sensitivity_grid.py` |
| **B** | Which region of the amplicon is a blind spot? | `detect_positional_bias.py` |
| **C** | What support threshold (`-e` or fractional) maximises F1? | `optimize_support_threshold.py` |
| **D** | How does Pindel compare to VarScan2, Mutect2, LoFreq on the same BAMs? | `compare_callers.py` |

**See [results/findings.md](results/findings.md) for the full results and interpretation.**

**Key finding:** LoFreq is the most reliable indel caller across both amplicons
(F1 0.73–0.74). Pindel requires per-amplicon insert-size calibration that cannot
be automated, and degrades severely on longer amplicons (LAGA_LC F1 0.35 vs 0.66
on tla). Recommend coupling LoFreq into the production amplicon pipeline.

## Status (August 2026)

All experiments complete on both datasets.

| Experiment | `tla` | `LAGA_LC` |
|---|---|---|
| A — sensitivity grid | ✓ | ✓ |
| B — positional blind spots | ✓ | ✓ |
| C — threshold sweep | ✓ | ✓ |
| D — caller comparison | ✓ | ✓ |

## Datasets

| ID | Description | Length | Reference |
|---|---|---|---|
| `tla` | Synthetic window of TLA D08 amplicon | 400 bp | `D08_HI_HS_AvB1.fa` pos 3000–3400 |
| `LAGA_LC` | Real LC amplicon, SAD0309 run | 711 bp | `SAD0309_P01_G04_LAGA-LC.fa` |

## Inputs

- Amplicon FASTA references (not tracked; see datasets above)
- ART MiSeq v3 error model (`MSv3`), paired-end, read lengths 75/150/250 bp
- No real FASTQ files — reads are simulated at controlled AF and depth

## Outputs

Result CSVs are in `results/expt{A,B,C,D}/<dataset>/` (tracked in git).
`nearmiss.csv` (~30–56 MB) is excluded from git.

```
results/
  exptA/<dataset>/   confusion_metrics.csv, detection_results.csv,
                     manifest.csv, sensitivity_by_length_af.csv
  exptB/<dataset>/   confusion_metrics.csv, detection_results.csv,
                     manifest.csv, position_sensitivity.csv
  exptC/<dataset>/   threshold_sweep.csv, af_gap.csv
  exptD/<dataset>/   caller_sensitivity.csv, caller_calls.csv
```

Run intermediates (BAMs, VCFs, mpileup, ~140 GB) are written to `runs/` which
is gitignored.

## Scripts

```
├── config.py                      shared constants and dataset geometry
├── detect_sensitivity_grid.py     Experiment A
├── detect_positional_bias.py      Experiment B
├── optimize_support_threshold.py  Experiment C
├── compare_callers.py             Experiment D
├── test_benchmark.py              guard tests
└── utils/
    ├── simulate.py                refs, ART simulation, AF mixing, run table
    ├── sim_pindel.py              trimmed Pindel runner (base -e 1 VCFs only)
    ├── align_call.py              bwa → samtools → sim_pindel → coverage plot
    ├── coverage_plot.py           per-condition depth PNG
    └── evaluate.py                match calls to truth, write result CSVs
```

## How to reproduce

```bash
conda activate indel_sim  # art_illumina, bwa, samtools, pindel, varscan, gatk4, lofreq, bcftools

# A/B accept -d {tla,LAGA_LC}; omit for both. D also accepts --callers, --force.

python detect_sensitivity_grid.py           # Exp A — run before C and D
python detect_positional_bias.py            # Exp B
python optimize_support_threshold.py        # Exp C (re-scores A's VCFs, no new sim)
python compare_callers.py                   # Exp D (re-calls A's BAMs, no new sim)
```

Conditions are resume-safe (existing `pindel-output/` or cached `.calls.tsv`
are skipped unless `--force`).

After running, copy result CSVs to the tracked `results/` folder:

```bash
BENCH=runs/sim_benchmark
for expt in A B C D; do
  for ds in tla LAGA_LC; do
    src="$BENCH/expt$expt/$ds/results"
    dst="results/expt$expt/$ds"
    [ -d "$src" ] && mkdir -p "$dst" && \
      find "$src" -maxdepth 1 -name '*.csv' ! -name 'nearmiss.csv' \
        -exec cp {} "$dst/" \;
  done
done
```

## Scoring

Matching is **exact** (`POS_TOLERANCE = SIZE_TOLERANCE = 0`): correct SVTYPE,
exact breakpoint position, exact size, support ≥ threshold. Right-type calls
within 15 bp of truth that fail exact matching go to `nearmiss.csv`.
AF=0 WT-only controls supply the FP baseline for precision.
