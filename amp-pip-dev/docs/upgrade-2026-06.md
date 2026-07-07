# Pipeline Upgrade — June 2026 (technical notes)

Technical record of the conda toolchain upgrade, the new structural-variant report,
and the bug fixes / tuning that followed. Each change carries its own **Tested** note
(the env upgrade is validated as output-neutral). Includes follow-up items for future
improvement. Audience: developers / maintainers.

**Dataset:** every experiment and update in this document was developed and validated
on the **2026-0609_TLA** run (3 samples), full plasmid reference `D08_HI_HS_AvB1.fa`.
Where old vs new environments are compared, both were run on this same dataset and
reference, so the only deliberate variable is the tool versions.

---

## 1. What changed

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

Reporting libs (pandas, numpy, openpyxl, plotly), `pindel=0.2.5b9`, and  `python-kaleido=0.2.1`
(chromium-bundled, for coverage-plot export) are unchanged in role.

**Tested:** Yes — and **output-neutral at the reporting level.** Old-env vs new-env
on the same reference gave a **byte-identical** `Results/snp.txt`, VarScan `*.snp.vcf`,
and Pindel `Summary.100.vcf` (all 3 samples). `*.indel.vcf`/`*.consensus.vcf` differ
only by ±1–2 read counts at a few positions with **no variant gained or lost**; the
only ≥2% change is the existing `6218 -GGCTCT` deletion shifting 0.03% (e.g. Lib1
30.85→30.82%), still called in both.

### B. Code: Pindel structural-variant report ([pipeline_combined.py](pipeline_combined.py))

Additive only — the `indel.txt` feature reads existing Pindel output; the tuned
Pindel invariants (insert size, `-E`, pindel2vcf sweep) are untouched. (The separate
`pindel_wrapper.py` thread fix is in §C.)

| Change | Detail |
|--------|--------|
| New helper `parse_pindel_summary()` | Parses each sample's existing `pindel-output/pindel-output.Summary.100.vcf` into a tidy DataFrame (Chrom, Position, Ref, Alt, SVTYPE, SVLEN, HOMLEN, GT, RefReads, AltReads, VarFreq). |
| New run-level report `Results/indel.txt` | Merges all samples' Pindel calls, laid out like `snp.txt`: grouped by clone (with replicate sample IDs) when clones are present, else by sample. Empty groups print "No structural variants reported by Pindel". |

Notes:
- Reads only existing per-sample Summary files; does not re-run or alter Pindel.
- Unlike `snp.txt` (which filters VarFreq ≥ 2%), `indel.txt` lists all calls from
  `Summary.100.vcf` and surfaces the AD read counts so the reviewer can judge each
  call (Pindel's GT field is not reliable at amplicon depth).

**Tested:** Yes. Validated parser + writer against existing Summary files (correct
grouped output) and syntax-checked in the new env.

### C. Thread settings (vs. the original `amp-pip`)

The per-stage thread counts in the wrappers were raised to 16 to match the 24-core
host (the original `amp-pip` used mixed, lower values). These are performance/
stability settings only — they do not change which variants are called.

| Wrapper | Original (`amp-pip`) | Now (`amp-pip-dev`) | Effect |
|---------|----------------------|----------------------|--------|
| `bwa_wrapper.py` | `bwa mem -t 6`, piped `samtools sort -@ 4` | `-t 16`, `-@ 16` | `bwa -t` causes ±1–2 read-count jitter (batch boundaries); sort thread count has none. |
| `samtools_wrapper.py` | `sort -@ 6`, `index -@ 6` | `-@ 16`, `-@ 16` | None (speed only; output deterministic). |
| `pindel_wrapper.py` | `pindel -T 100` | `-T 16` | **Fix:** `-T 100` oversubscribed a 24-core host (thread thrashing / potential stalls). Output unchanged. |
| `varscan_wrapper.py` | — | — | No thread flag; unchanged. |

The only result-relevant one is `bwa mem -t`, and only as sub-threshold read-count
jitter: comparing two runs that differ in thread count only (same env and reference),
BWA-MEM batch boundaries shift with `-t` so per-position read counts move by ±1–2
while the called variants stay identical.

**Tested:** Yes — the ±1–2 jitter above stays below every reporting threshold; no
variant gained or lost.

### D. Bug fixes — fastp report paths & coverage-plot span

Correctness fixes found while validating on real MiSeq (`_R1_001.fastq.gz`) naming.
None change which variants are called.

| Change | File | Detail |
|--------|------|--------|
| fastp report path no longer clobbers input | `fastp_wrapper.py` | The old `input.replace(".R1.fastq.gz", ".fastp.html")` **silently no-ops** on `_R1_001.fastq.gz`-style names, so fastp wrote its HTML/JSON report over the raw input FASTQ. Now the FASTQ extension (`.fastq.gz`/`.fq.gz`/`.fastq`/`.fq`) is stripped before appending `.fastp.html`/`.fastp.json`. |
| Coverage plot spans the full reference | `plot_coverage_plots.py` | The plot was drawn only from positions present in the consensus, so it started partway in and stopped short. New helper `get_reference_length()` reads the reference FASTA; `read_consensus_file()` now pads zero-coverage positions before the first and after the last covered base up to the reference length. |
| Per-sample coverage-plot filename | `plot_coverage_plots.py` | Output renamed from a fixed `coverage_plot.jpeg` to `<cell_line_id>_coverage_plot.jpeg` so multi-sample runs no longer overwrite one image. |
| Pass reference through to the plotter | `pipeline.py` | `plot_coverage_plots_manager()` now receives `reference_file_path` so it can compute the reference length for the padding above. |

### E. mpileup depth cap removed — `samtools mpileup -d 0`

`samtools_wrapper.py` now passes an explicit `-d 0` (uncapped max depth) instead
of relying on the samtools default `-d 8000`. **This deliberately changes output**
(unlike §C/§D): the 8000 cap discards ~98% of reads on this deep amplicon and
inflates a handful of borderline positions over the 2% cutoff — i.e. VarScan SNP
false positives. `-d 0` feeds VarScan the true read set and removes those spurious
≥2% calls, making the SNP report more stringent.

Scope: correct for this small (~10 kb) ultra-deep amplicon assay; **do not** copy
`-d 0` to WGS/large-panel pipelines (OOM risk at high-coverage loci). Full
investigation and per-library numbers: `experiments/2026-06-mpileup-depth-cap/`.

---

## 2. Future improvements / follow-ups

1. **Thread-independent alignment (optional).** BWA-MEM read-count jitter comes from
   batch boundaries shifting with thread count. `bwa mem -K 10000000` makes batching
   deterministic and removes the ±1–2 read differences seen in the consensus VCFs.
   Cosmetic only — does not affect reported variants.
2. **Consider a VarFreq display threshold for `indel.txt`.** Currently it shows all
   `Summary.100.vcf` calls (including low-AD ones). If reviewers prefer parity with
   `snp.txt`, add an optional frequency filter — but keep raw AD visible, since
   Pindel GT is unreliable at amplicon depth.
3. **Migrate `pileup2*` → `mpileup2*` (needs an experiment first).** We now generate
   the `.mpileup` file but still call VarScan's **legacy** `pileup2snp`/`pileup2indel`/
   `pileup2cns`. Whether to lift to the modern `mpileup2*` callers should be decided by
   a follow-up experiment (extend `experiments/2026-06-varscan-interface/`), which found
   the two interfaces are **not** equivalent. **Caution:** the interfaces have different
   default reporting `--min-var-freq` (pileup2snp 0.01 vs mpileup2snp 0.2) and emit
   different output formats — any migration must force matched params and update the
   downstream CONS-format parsers.
4. **Reorganize the run output structure for transfer/storage.** Large intermediate
   files (unsorted/sorted `.bam`, `.mpileup`) are never manually inspected but dominate
   run size and slow user file transfer. Consider a cleaner layout that separates
   reviewer-facing reports (`snp.txt`, `indel.txt`, coverage plots) from bulky
   intermediates, so the latter can be pruned or archived without hunting through
   per-sample folders.
5. **Tune the hard-coded Pindel insert size (`398`).** `pindel_wrapper.py` hard-codes
   `398` as the insert-size config. This should be experimented against the value
   implied by the actual read length / library insert size for the most accurate
   indel calls, and the chosen value documented.
---

## 3. Files touched

| File | Change | Tested |
|------|--------|--------|
| [envs/variant-calling.yml](../envs/variant-calling.yml) | New env definition (Python 3.10 + upgraded tools) | Yes (env built + full run) |
| [pipeline_combined.py](../pipeline_combined.py) | Added Pindel → `indel.txt` report; renamed from `pindel.txt`; updated docstrings | Yes |
| [bwa_wrapper.py](../bwa_wrapper.py) | `bwa mem -t 6→16`, piped `samtools sort -@ 4→16` | Yes |
| [samtools_wrapper.py](../samtools_wrapper.py) | `sort`/`index` `-@ 6→16`; **`mpileup -d 0`** (uncapped, §1.E) | Yes |
| [pindel_wrapper.py](../pindel_wrapper.py) | `-T 100→16` (thread fix); tuned invariants unchanged | Yes |
| [fastp_wrapper.py](../fastp_wrapper.py) | Report path no longer clobbers input FASTQ on `_R1_001` names (§1.D) | Yes |
| [plot_coverage_plots.py](../plot_coverage_plots.py) | Pad coverage to full reference length; per-sample plot filename; new `get_reference_length()` (§1.D) | Yes |
| [pipeline.py](../pipeline.py) | Pass `reference_file_path` to `plot_coverage_plots_manager()` (§1.D) | Yes |

## 4. Reproduce a run

```bash
cd ngs/amp-pip-dev
conda run --no-capture-output -n variant-calling \
  python pipeline_combined.py /abs/path/to/RUN_DIR
```
Run from the repo root so the `../amino-anno/from_varscan.py` sibling path resolves.
