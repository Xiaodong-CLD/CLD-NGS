# Cell Line Development amplicon-seq pipeline

Amplicon-seq variant calling for cell-line development. Takes paired-end FASTQ
files through trimming, alignment, and variant calling, then produces per-sample
reports plus replicate-aware run-level summaries of SNPs and indels.

**Pipeline stages:** trim (fastp) → align (BWA-MEM) → pileup (samtools) →
SNV calling (VarScan) + indel calling (Pindel, LoFreq, VarScan) → per-sample
reports and a replicate-aware run-level summary.

> This pipeline runs **locally** (not on a Slurm HPC). Ignore any older Slurm /
> `job.script` instructions.

---

## Dependencies

Install the conda environment once before running:

```bash
conda env create -f envs/variant-calling.yml
conda activate variant-calling
```

Key tools pinned in `envs/variant-calling.yml`:

| Tool | Version |
|------|---------|
| Python | 3.10.20 |
| BWA | 0.7.19 |
| samtools | 1.23.1 |
| fastp | 1.3.5 |
| VarScan | 2.4.6 |
| Pindel | 0.2.5b9 |
| LoFreq | 2.1.5 |
| OpenJDK | 17.0.18 |
| pandas | 2.3.3 |
| openpyxl | 3.1.5 |
| plotly | 6.8.0 |

This pipeline also calls `../amino-anno/from_varscan.py` (sibling repository) for
amino-acid annotation of SNP calls. That repo must be present at the same directory
level as `amp-pip-dev`.

---

## Sample sheet format

The configuration spreadsheet has one row per sample. The first column is used as
a unique sample identifier (index). Required columns — exact spelling matters:

| Column | Required | Meaning |
|--------|----------|---------|
| `fastq (forward)` | Yes | Forward (R1) FASTQ filename inside `--fastq-dir` |
| `fastq (reverse)` | Yes | Reverse (R2) FASTQ filename inside `--fastq-dir` |
| `GOI reference file` | Yes | Reference FASTA **basename** — must match the file passed as `--reference` |
| `Clone ID (optional)` | No | If present and filled, groups replicates under a clone for the run-level summary. The column may be omitted entirely or left blank to report per-sample. |

> **Important:** `GOI reference file` must contain the exact basename of the
> `--reference` file. If you pass `--reference /path/to/goi.fa`, every row must
> contain `goi.fa` in that column.

---

## Running

```bash
conda activate variant-calling
cd /path/to/amp-pip-dev

python run_pipeline.py \
    --sample-sheet  /path/to/sample_sheet.xlsx \
    --fastq-dir     /path/to/Fastq \
    --reference     /path/to/goi.fa \
    --output-dir    /path/to/output
```

Windows (`Z:`) paths are accepted when the Z: drive is mounted at `/mnt/z`.

**Alternatively**, if you have already built the legacy input directory layout
(`Fastq/`, `reference/`, `sample_xlsx/` inside a single folder), call
`pipeline_combined.py` directly:

```bash
python pipeline_combined.py <input_directory>
```

---

## Outputs

All outputs are written under `<output-dir>/`:

```
<output-dir>/
├─ Results/
│  ├─ snp.txt                      run-level SNP report (VarFreq ≥ 1%)
│  ├─ indel.txt                    run-level indel report (≥ 1% AF, all callers)
│  └─ <sample>/
│     ├─ *.fastp.html / .json      fastp QC reports
│     ├─ coverage-plot-output/     per-sample coverage plot
│     ├─ varscan-output/           SNP VCF, indel VCF, amino-annotated SNP VCF
│     ├─ pindel-output/            pindel-output.Summary.vcf
│     └─ lofreq-output/            <sample>.lofreq.indel.vcf
└─ tmp/
   └─ <sample>/                    intermediates (BAM, mpileup, raw VCFs)
                                   safe to delete after reviewing Results/
```

- **`Results/snp.txt`** — run-level SNP summary (VarScan, variants ≥ 1% AF).
  Replicate-aware when `Clone ID (optional)` is filled for every sample: a
  variant is reported only when it appears in ≥ 1 replicate of a clone. Falls
  back to per-sample reporting when clone IDs are absent.
- **`Results/indel.txt`** — run-level multi-caller indel summary (Pindel +
  VarScan + LoFreq). All positions with ≥ 1% read support from any caller are
  reported, with a `source` column indicating which caller(s) detected them.
  Grouped by clone (or by sample), same layout as `snp.txt`.
- **`Results/<sample>/varscan-output/*.snp.amino_annotated.vcf`** — VarScan SNP
  calls with amino-acid consequence annotation from `amino-anno`.

---

## Thresholds

| Filter | Value | Location |
|--------|-------|----------|
| SNP minimum coverage | 1500 reads | `varscan_wrapper.py` |
| SNP report frequency | ≥ 1% VarFreq | `pipeline_combined.py` |
| Indel report frequency | ≥ 1% AF (any caller) | `pipeline_combined.py` |
| samtools pileup depth cap | uncapped (`-d 0`) | `samtools_wrapper.py` |

The uncapped pileup depth is intentional for deep amplicon sequencing; do not
copy this setting to WGS or large-panel pipelines.

---

## Indel calling — known limitations

**Pindel insert size is amplicon-specific.** The insert size passed to Pindel
(`-i`, hardcoded as `398` in `pindel_wrapper.py`) is calibrated for ~400 bp
TLA-geometry amplicons. It yields **zero true calls** on longer amplicons such
as LAGA-LC (711 bp), where the working range is ~155–175. If you are running
on a new amplicon geometry and Pindel returns no indels, check whether `-i 398`
is appropriate for your fragment size. See
`expt_summer2026/amp-pip-dev/2026-07_pindel_indel_sensitivity/results/findings.md`
for the calibration results.

**LoFreq is the most reliable indel caller.** Benchmarking on both the TLA and
LAGA-LC amplicons showed LoFreq achieves F1 0.73–0.74 consistently, compared to
Pindel's 0.35–0.66 (geometry-dependent) and VarScan's 0.40–0.54. LoFreq calls
should be the primary source for indel interpretation when Pindel and VarScan
disagree or return nothing.

---

## What's new (August 2026)

- **`run_pipeline.py`** — new user-facing launcher; replaces the older
  `pipeline_runner.py` (now in `archive/`).
- **LoFreq indel caller** added alongside Pindel and VarScan.
- **`Results/indel.txt`** — new multi-caller indel report.
- **`Clone ID (optional)` column** is now truly optional; the column may be
  omitted from the sample sheet entirely.
- **Streamlit launcher** — pipeline is available through `pipeline_app`.

See [docs/upgrade-2026.md](docs/upgrade-2026.md) for the full technical
change log (June and August 2026).

---

## Repository layout

```
amp-pip-dev/
├─ run_pipeline.py          user-facing launcher (start here)
├─ pipeline_combined.py     primary orchestrator; reads xlsx, calls pipeline.py
├─ pipeline.py              single-sample runner
├─ bwa_wrapper.py           BWA-MEM alignment
├─ samtools_wrapper.py      samtools sort / pileup
├─ fastp_wrapper.py         adapter trimming and QC
├─ varscan_wrapper.py       VarScan SNV and indel calling
├─ pindel_wrapper.py        Pindel structural-variant calling + pindel2vcf
├─ lofreq_wrapper.py        LoFreq indel calling
├─ filter_deletions.py      filters RPL/HOMLEN entries from Pindel _D output
├─ merge_fastq.py           FASTQ pair discovery helper
├─ plot_coverage_plots.py   per-sample coverage plot generation
├─ envs/                    conda environment files
├─ docs/                    upgrade notes and user manual
│  ├─ upgrade-2026.md       technical change log
│  └─ User Manual Amp Var Calling.pdf
└─ archive/                 scripts no longer used by the pipeline
   ├─ pipeline_runner.py    old multi-sample launcher (superseded by run_pipeline.py)
   ├─ all_2_fasta.py        generic FASTA converter (unrelated utility)
   ├─ dependency_simplify.py  conda env utility
   └─ manual_tools/         preflight and ad-hoc tools (not part of the pipeline run)
```

---

## Further reading

- [docs/upgrade-2026.md](docs/upgrade-2026.md) — full technical change log
- [docs/User Manual Amp Var Calling.pdf](docs/User%20Manual%20Amp%20Var%20Calling.pdf) — end-user manual
- `expt_summer2026/amp-pip-dev/` — benchmarking experiments (mpileup depth, Pindel insert size, caller comparison)
