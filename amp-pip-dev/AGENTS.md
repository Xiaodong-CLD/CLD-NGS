# AGENTS.md — amp-pip (amplicon-seq variant calling)

> Status: **under active development.** Behavior reflects the current code, which sometimes diverges from [README.md](README.md). When they disagree, trust the code and flag the mismatch.
> Assumption: the pipeline now runs **locally, not on a Slurm-based HPC**. Ignore Slurm/`job.script` instructions in the README.

## 1. Purpose
Amplicon-seq variant calling for cell-line development: paired-end FASTQ → trim → BWA-MEM align → samtools pileup → VarScan SNV/indel + Pindel structural variants → per-sample reports and a replicate-aware SNP summary.

## 2. Entry points
- [pipeline.py](pipeline.py) — single sample: `pipeline.py -r REF -f1 R1 -f2 R2 [-nt]`.
- [pipeline_combined.py](pipeline_combined.py) — **primary runner:** `pipeline_combined.py INPUT_DIR`. Reads the xlsx config, runs `pipeline.py` per sample, adds amino-acid annotation, writes the summary.
- [pipeline_runner.py](pipeline_runner.py) — runs `pipeline.py` over a directory of sample subfolders.
- [Run_pindel_only.py](Run_pindel_only.py) — Pindel-only re-run. [check_pipeline_tools_requirement.py](check_pipeline_tools_requirement.py) — preflight tool check.

## 3. Inputs / outputs (per `pipeline_combined.py`)
- Inputs under `INPUT_DIR/`: `sample_xlsx/*.xlsx` (first match), `Fastq/`, `reference/`.
- xlsx columns (exact strings): `fastq (forward)`, `fastq (reverse)`, `GOI reference file`, `Clone ID (optional)`.
- Outputs: `INPUT_DIR/Results/<sample>/` (incl. `varscan-output/`, `pindel-output/`, `coverage-plot-output/`); summary `INPUT_DIR/Results/snp.txt`.
- Shells out to `../amino-anno/from_varscan.py` (sibling repo) — keep that path intact.

## 4. Invariants — do not change silently
- VarScan SNP `--min-coverage 1500` ([varscan_wrapper.py](varscan_wrapper.py)) and summary `VarFreq >= 2%` filter ([pipeline_combined.py](pipeline_combined.py)) — biological thresholds.
- Pindel ([pindel_wrapper.py](pindel_wrapper.py)): insert size `398`, `-E 0.9`, pindel2vcf sweep `[15, 30, 50, 75, 100]`; summary built only from `*.100.vcf`.
- Sample naming: `name()` splits the forward FASTQ on `_L`; renamed reads use `.R1`/`.R2`.
- Filename suffixes (`.bam`, `.sorted.bam`, `.mpileup`, `.trimmed.fastq.gz`, `.snp.vcf`, `.consensus.vcf`) are an implicit contract between stages.
- Resume-by-existence: wrappers skip work if output exists; delete partial files to force a re-run.

## 5. Rules for code changes
- Make the smallest useful change; don't refactor unrelated stages.
- Prefer simple Python helpers over new classes or abstractions.

## 6. Do not refactor (without explicit ask)
- The `subprocess` + `shell=True` command strings — they encode tuned tool flags.
- The filename `.replace()`/`.split()` conventions and skip-if-exists logic (see §4).
- The `../amino-anno/from_varscan.py` dependency and `sample_xlsx/Fastq/reference/Results` layout.
- Don't edit code to match the README's paths — confirm intent first.
