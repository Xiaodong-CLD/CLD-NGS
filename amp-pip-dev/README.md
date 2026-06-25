# Cell Line Development amplicon-seq pipeline

Amplicon-seq variant calling for cell-line development. Takes paired-end FASTQ
files through trimming, alignment, and variant calling, then produces per-sample
reports plus replicate-aware run-level summaries of SNPs and structural variants.

**Pipeline stages:** trim (fastp) → align (BWA-MEM) → pileup (samtools) →
SNV/indel calling (VarScan) + structural variants (Pindel) → per-sample reports
and a replicate-aware summary.

> This pipeline runs **locally** (not on a Slurm HPC). Ignore any older Slurm /
> `job.script` instructions.

---

## What's new (June 2026)

- **Updated tool environment** (`variant-calling`): Python 3.10 with refreshed
  BWA, samtools, fastp, and VarScan. Variant-calling results are unchanged versus
  the previous environment (validated — see [docs/upgrade-2026-06.md](docs/upgrade-2026-06.md)).
- **New structural-variant report `Results/indel.txt`**: a run-level summary of
  Pindel insertion/deletion calls, grouped by clone (or by sample when no Clone ID
  is given) — mirroring the layout of `snp.txt`.

See [docs/upgrade-2026-06.md](docs/upgrade-2026-06.md) for the full technical write-up.

---

## Setup

Create the conda environment from the provided file (one time):

```bash
conda env create -f envs/variant-calling.yml
```

Activate it before running:

```bash
conda activate variant-calling
```

---

## Preparing a run

Each run lives in its own input directory with this layout:

```
<input directory>/
├─ sample_xlsx/      # one .xlsx configuration spreadsheet (first match is used)
├─ Fastq/            # all forward and reverse FASTQ files referenced in the sheet
└─ reference/        # all GOI reference .fasta files referenced in the sheet
```

The configuration spreadsheet has one row per sample, with these exact columns:

| Column | Meaning |
|--------|---------|
| `fastq (forward)` | Forward (R1) FASTQ filename in `Fastq/` |
| `fastq (reverse)` | Reverse (R2) FASTQ filename in `Fastq/` |
| `GOI reference file` | Reference `.fasta` filename in `reference/` |
| `Clone ID (optional)` | Groups replicates under a clone; optional |

---

## Running

Run from the repository root (so the sibling `../amino-anno/from_varscan.py`
amino-acid annotation step resolves):

```bash
python pipeline_combined.py <input directory>
```

---

## Outputs

All outputs are written under `<input directory>/Results/`:

- **`Results/<sample>/`** — per-sample outputs (`varscan-output/`, `pindel-output/`,
  `coverage-plot-output/`). The VarScan SNP report also gets an amino-acid
  annotation file.
- **`Results/snp.txt`** — run-level SNP summary (variants ≥ 2% frequency). If the
  `Clone ID` column is filled for every sample, the summary is replicate-aware
  (a variant is reported only when it appears in more than one replicate of a
  clone). If any sample lacks a Clone ID, the summary is split out by sample.
- **`Results/indel.txt`** — run-level Pindel structural-variant summary, grouped
  by clone (or by sample), in the same layout as `snp.txt`.
