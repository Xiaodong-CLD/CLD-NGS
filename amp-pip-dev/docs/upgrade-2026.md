# Pipeline Upgrade — 2026 (technical notes)

Technical record of all changes made to `amp-pip-dev` in 2026. Audience:
developers / maintainers.

**Dataset for June 2026 validation:** all experiments and comparisons in the
June section were run on the **2026-0609_TLA** dataset (3 samples), full plasmid
reference `D08_HI_HS_AvB1.fa`. Old-env vs new-env comparisons use this as the
only varying factor.

---

## June 2026: toolchain upgrade, initial indel report, bug fixes

### A. Conda environment toolchain upgrade

New environment `variant-calling` (from [envs/variant-calling.yml](envs/variant-calling.yml))
replaces the previous `variant-calling` environment.

| Tool      | Old (baseline) | New (upgraded) | Reason |
|-----------|----------------|----------------|--------|
| python    | 3.8.17         | 3.10.20        | Standardize on Python 3.10 (standing preference); 3.8 is end-of-life. |
| bwa       | 0.7.17         | 0.7.19         | Latest 0.7.x bug-fix release. |
| samtools  | 1.17           | 1.23.1         | Current stable; htslib improvements. |
| fastp     | 0.23.4         | 1.3.5          | Major version jump (0.x → 1.x); current stable trimmer. **Largest behavior-change risk** — see validation. |
| varscan   | 2.4.4          | 2.4.6          | Latest 2.4.x. |
| openjdk   | (implicit)     | 17             | Pinned modern JRE for VarScan (Java). |

Reporting libs (pandas, numpy, openpyxl, plotly), `pindel=0.2.5b9`, and
`python-kaleido=0.2.1` (chromium-bundled, for coverage-plot export) are
unchanged in role.

**Tested:** Yes — and **output-neutral at the reporting level.** Old-env vs new-env
on the same reference gave a **byte-identical** `Results/snp.txt`, VarScan `*.snp.vcf`,
and Pindel `Summary.100.vcf` (all 3 samples). `*.indel.vcf`/`*.consensus.vcf` differ
only by ±1–2 read counts at a few positions with **no variant gained or lost**.

### B. Pindel structural-variant report (`pipeline_combined.py`)

Additive only — reads existing Pindel output; tuned Pindel invariants (insert
size, `-E`, pindel2vcf sweep) are untouched in this release.

| Change | Detail |
|--------|--------|
| New helper `parse_pindel_summary()` | Parses each sample's `pindel-output/pindel-output.Summary.100.vcf` into a tidy DataFrame (Chrom, Position, Ref, Alt, SVTYPE, SVLEN, HOMLEN, GT, RefReads, AltReads, VarFreq). |
| New run-level report `Results/indel.txt` | Merges all samples' Pindel calls, laid out like `snp.txt`: grouped by clone when clones are present, else by sample. Empty groups print "No structural variants reported by Pindel". |

**Tested:** Yes. Validated parser + writer against existing Summary files and
syntax-checked in the new env.

### C. Thread settings

Per-stage thread counts raised to 16 to match the 24-core host.

| Wrapper | Original (`amp-pip`) | Now (`amp-pip-dev`) | Effect |
|---------|----------------------|----------------------|--------|
| `bwa_wrapper.py` | `bwa mem -t 6`, piped `samtools sort -@ 4` | `-t 16`, `-@ 16` | `bwa -t` causes ±1–2 read-count jitter (batch boundaries); sort thread count has none. |
| `samtools_wrapper.py` | `sort -@ 6`, `index -@ 6` | `-@ 16`, `-@ 16` | None (speed only). |
| `pindel_wrapper.py` | `pindel -T 100` | `-T 16` | **Fix:** `-T 100` oversubscribed a 24-core host. Output unchanged. |

### D. Bug fixes — fastp report paths & coverage-plot span

| Change | File | Detail |
|--------|------|--------|
| fastp report path no longer clobbers input | `fastp_wrapper.py` | Old `input.replace(".R1.fastq.gz", ".fastp.html")` silently no-ops on `_R1_001.fastq.gz`-style names, causing fastp to overwrite the raw FASTQ. Fixed: strip the FASTQ extension before appending `.fastp.html`/`.fastp.json`. |
| Coverage plot spans the full reference | `plot_coverage_plots.py` | Plot started partway in and stopped short. New `get_reference_length()` reads the reference FASTA; `read_consensus_file()` pads zero-coverage positions at both ends up to reference length. |
| Per-sample coverage-plot filename | `plot_coverage_plots.py` | Output renamed from `coverage_plot.jpeg` to `<cell_line_id>_coverage_plot.jpeg` so multi-sample runs no longer overwrite one image. |
| Pass reference through to the plotter | `pipeline.py` | `plot_coverage_plots_manager()` now receives `reference_file_path`. |

### E. mpileup depth cap removed — `samtools mpileup -d 0`

`samtools_wrapper.py` now passes `-d 0` (uncapped depth) instead of the samtools
default `-d 8000`. The 8000 cap discards ~98% of reads on deep amplicons and
inflates a handful of borderline positions over the 2% cutoff (VarScan SNP false
positives). `-d 0` feeds VarScan the true read set and removes those spurious calls.

Scope: correct for this small (~10 kb) ultra-deep amplicon assay; **do not** copy
`-d 0` to WGS/large-panel pipelines (OOM risk at high-coverage loci). Full
investigation: `expt_summer2026/amp-pip-dev/2026-06_amp-pip-update/2026-06-mpileup-depth-cap/`.

### Files touched (June 2026)

| File | Change |
|------|--------|
| `envs/variant-calling.yml` | Toolchain upgrade (Python 3.10, fastp 1.x, samtools 1.23.1, etc.) |
| `bwa_wrapper.py` | Thread count raised to 16 |
| `samtools_wrapper.py` | Thread count raised to 16; `-d 0` added to mpileup |
| `pindel_wrapper.py` | `-T 100` → `-T 16` |
| `fastp_wrapper.py` | Report path fix |
| `plot_coverage_plots.py` | Reference-length padding; per-sample filename |
| `pipeline.py` | Pass `reference_file_path` to plotter |
| `pipeline_combined.py` | Added `parse_pindel_summary()`, initial `Results/indel.txt` writer |

---

## August 2026: LoFreq, multi-caller indels, new launcher

### A. User-facing entry point (`run_pipeline.py`)

Added `run_pipeline.py` as the new user-facing launcher. Previously users had
to manually build the expected input directory structure before calling
`pipeline_combined.py`.

**New interface:**

```bash
python run_pipeline.py \
    --sample-sheet  samples.xlsx \
    --fastq-dir     /path/to/Fastq \
    --reference     goi.fa \
    --output-dir    /path/to/output
```

**How it works (compatibility shim):**

1. Calls `tools/prepare_pipeline_input.sh` to create symlinks inside
   `<output-dir>`:

   ```
   <output-dir>/Fastq/          ← symlinks to source FASTQ files
   <output-dir>/sample_xlsx/    ← symlink to sample sheet
   <output-dir>/reference/      ← symlink to reference FASTA
   ```

2. Calls `pipeline_combined.py <output-dir>` using the prepared directory as
   both input source and output destination.

No FASTQ files are copied. The internal pipeline path logic is **unchanged**.

### B. `tools/prepare_pipeline_input.sh`

Copied from `expt_summer2026/amp-pip-dev/scripts/prepare_pipeline_input.sh`
into the pipeline `tools/` directory so `run_pipeline.py` can locate it
relative to itself.

### C. Sample order

Sample order is already preserved. `pd.read_excel(..., index_col=0).to_dict(orient='index')`
returns keys in spreadsheet row order (Python 3.7+ dict insertion order), so
`Results/snp.txt` and `Results/indel.txt` iterate samples in the order they
appear in the sample sheet. No change was needed.

### D. Tests

Added `tests/test_run_pipeline.py` with two tests:

- `test_run_pipeline_calls_prepare_sh_with_correct_args` — asserts that
  `run_pipeline.py` calls the bash script and `pipeline_combined.py` with the
  correct positional arguments.

- `test_run_pipeline_produces_same_results_as_direct_call` — end-to-end smoke
  test confirming that the compatibility layer produces the same `Results/` layout
  as calling `pipeline_combined.py` directly.

### E. LoFreq indel caller (`lofreq_wrapper.py`, `pipeline.py`, `pipeline_combined.py`)

Added LoFreq as a third indel caller alongside Pindel and VarScan.

| Change | Detail |
|--------|--------|
| New `lofreq_wrapper.py` | Runs `lofreq call-parallel` on the sorted BAM with `--call-indels`. Outputs `<sample>.lofreq.indel.vcf` under `tmp/<sample>/lofreq-output/`. |
| `pipeline.py` — `call_indels()` | Now calls both `pindel_runner()` and `lofreq_runner()`. |
| `pipeline_combined.py` — per-sample loop | Calls `parse_lofreq_indels()` on the LoFreq VCF, appends to per-sample lists, copies final VCF to `Results/<sample>/lofreq-output/`. |

**Tested:** LoFreq runs successfully on real MiSeq amplicon data and produces
indel VCFs in the expected format.

### F. Multi-caller indel aggregation (`pipeline_combined.py`)

`Results/indel.txt` was Pindel-only in June 2026. It now merges calls from all
three callers.

| Change | Detail |
|--------|--------|
| New helpers `parse_varscan_indels()`, `parse_lofreq_indels()` | Parse VarScan indel TSV and LoFreq indel VCF respectively into tidy DataFrames with consistent columns. |
| New helper `merge_indel_callers()` | Merges on `(SampleID, Chrom, Position)`; all positions detected by any caller at ≥ 1% AF are included. Adds a `source` column showing which caller(s) found it (e.g. `pindel`, `pindel,lofreq`). |
| Updated `Results/indel.txt` writer | Iterates the merged DataFrame grouped by clone/sample (same layout as `snp.txt`). NaN caller-freq cells replaced with `-` for readability. |
| Minimum support thresholds | Pindel: AD alt / (ref + alt) ≥ 1% (`PINDEL_MIN_AD_FRAC = 0.01`). VarScan: VarFreq ≥ 1%. LoFreq: AF ≥ 1%. |

**Rationale:** each caller's 1% threshold filters sequencing noise independently.
The `source` column lets the reviewer see at a glance whether a call is
corroborated by multiple callers or is single-caller only.

**Note on Pindel AF:** Pindel's `AD` field over-reports true allele fraction
approximately 4× at amplicon depth (sensitivity benchmarks:
`expt_summer2026/amp-pip-dev/2026-07_pindel_indel_sensitivity/`). The 1% `AD`
fraction threshold is therefore not equivalent to a true 1% AF cutoff and should
be interpreted with the `source` annotation — corroboration from VarScan or
LoFreq gives higher confidence.

### G. `pindel2vcf` — single call with in-Python filter (replaces `-e` sweep)

Previously `pindel_wrapper.py` ran `pindel2vcf` five times per SV type
(`-e 15/30/50/75/100`); only the `-e 100` output (`Summary.100.vcf`) was ever
used downstream. The `-e` flag is a post-hoc filter on the alt-allele supporting
read count (verified: the `-e 100` VCF is exactly the `-e 1` VCF with rows whose
`AD` alt < 100 dropped), so the intermediate VCFs were written and never read.

The pipeline now runs a single `pindel2vcf` call at `-e 1` and applies
`PINDEL_MIN_AD_FRAC = 0.01` as an in-Python fraction filter inside
`merge_indel_callers()`. This removes redundant disk writes and makes the
threshold explicit and adjustable without re-running Pindel.

### H. Streamlit launcher integration (`pipeline_app/pipelines/amp.py`)

The pipeline is now available through the `pipeline_app` Streamlit interface.
The form accepts the four user-facing inputs (`--sample-sheet`, `--fastq-dir`,
`--reference`, run name) and builds a `run_pipeline.py` command. Output is
written to `pipeline_outputs/amp-pip/<run_name>/`.

### Reference filename — important note

`run_pipeline.py` symlinks `--reference` into `<output-dir>/reference/` under
its original basename. `pipeline_combined.py` then constructs the per-sample
reference path as:

```
<output-dir>/reference/<GOI reference file column value>
```

**The `GOI reference file` column in the sample sheet must therefore contain the
exact basename of the file passed as `--reference`.** A mismatch causes BWA to
fail with "No such file or directory" even though fastp succeeds (fastp reads
FASTQs, not the reference). This is a known limitation of the symlink shim — see
the future-work item below.

### Tests

All tests (`tests/test_routing.py`, `tests/test_run_pipeline.py`) pass with the
August 2026 changes. The `test_run_pipeline.py` mock seeds the LoFreq output
directory; `test_routing.py` verifies that `lofreq-output/` is copied to
`Results/<sample>/`.

### Files touched (August 2026)

| File | Change |
|------|--------|
| `run_pipeline.py` | New user-facing entry point |
| `tools/prepare_pipeline_input.sh` | Copied from `expt_summer2026/amp-pip-dev/scripts/` |
| `lofreq_wrapper.py` | New — LoFreq caller wrapper |
| `pindel_wrapper.py` | `-e` sweep removed; single `pindel2vcf -e 1` call |
| `pipeline.py` | `call_indels()` now calls both Pindel and LoFreq |
| `pipeline_combined.py` | Added `parse_varscan_indels()`, `parse_lofreq_indels()`, `merge_indel_callers()`; updated indel report writer; LoFreq output copied to `Results/` |
| `pipeline_app/pipelines/amp.py` | Updated to 4-input form pointing at `run_pipeline.py` |
| `tests/test_run_pipeline.py` | New end-to-end smoke tests for the launcher |
| `README.md` | Updated input instructions, running instructions, output layout |
| `docs/upgrade-2026.md` | This file |

---

## Future improvements / open items

1. **Migrate `pileup2*` → `mpileup2*` (needs dedicated experiment, proceed with
   caution).** VarScan legacy `pileup2snp`/`pileup2indel`/`pileup2cns` are still
   used. The two interfaces are not equivalent: they differ in default
   `--min-var-freq` (pileup2snp 0.01 vs mpileup2snp 0.2) and emit different output
   formats. Any migration must force matched params and update the downstream
   CONS-format parsers. Investigate in
   `expt_summer2026/amp-pip-dev/2026-06_amp-pip-update/2026-06-varscan-interface/`
   before touching production code.

2. **Tune hard-coded Pindel insert size (`398`).** `pindel_wrapper.py` hard-codes
   `398`. Should be validated against the actual library insert size and documented.

3. **Reorganize run output for transfer/storage.** Large intermediates (sorted BAMs,
   `.mpileup`) dominate run size and slow file transfer. Consider separating
   reviewer-facing reports (`snp.txt`, `indel.txt`, coverage plots) from bulky
   intermediates so they can be pruned or archived independently.

4. **Remove the symlink shim — refactor `pipeline_combined.py` to accept paths
   directly.** The `run_pipeline.py` → `prepare_pipeline_input.sh` →
   `pipeline_combined.py` chain is a compatibility shim. The goal is to pass
   `--sample-sheet`, `--fastq-dir`, `--reference`, and `--output-dir` directly into
   `pipeline_combined.py` so that:
   - `run_pipeline.py` and `tools/prepare_pipeline_input.sh` can be removed.
   - The `Fastq/`, `sample_xlsx/`, and `reference/` symlink directories no longer
     appear inside `--output-dir`.
   - The reference filename mismatch issue is eliminated — the script resolves the
     reference path from `--reference` directly rather than reconstructing it from
     the sample sheet column.
   - Inputs and outputs live in separate, user-chosen locations.

   Affected code in `pipeline_combined.py`: the `config_f` glob, the
   `fq_1_path`/`fq_2_path` construction, and the `reference_path` construction
   blocks in the main script body.
