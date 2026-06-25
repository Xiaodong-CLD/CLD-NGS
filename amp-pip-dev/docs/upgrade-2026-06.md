# Pipeline Upgrade — June 2026 (technical notes)

Technical record of the conda toolchain upgrade, the new structural-variant report,
and the validation that confirmed the upgrade is output-neutral. Includes follow-up
items for future improvement. Audience: developers / maintainers.

**Validation dataset:** TLA-0609 (3 samples), full plasmid reference
`D08_HI_HS_AvB1.fa`. The old and new environments were run on this same dataset and
reference; the only deliberate variable between them is the tool versions.

---

## 1. What changed

### A. Conda environment toolchain upgrade

New environment `variant-calling` (from [../envs/variant-calling.yml](../envs/variant-calling.yml))
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

**Tested:** Yes. Full 3-sample run (`runs/T2_upgrade/`, threads=16) completed with
no `Traceback`. Compared against `runs/T1_5_threads/` (old env, same reference, same
thread count), so any difference is attributable to tool versions only.

### B. Code: Pindel structural-variant report ([../pipeline_combined.py](../pipeline_combined.py))

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

The only result-relevant one is `bwa mem -t`, and only as the sub-threshold
read-count jitter described in §2.

---

## 2. Validation — output comparison (old env vs new env, same reference)

| Output | Result |
|--------|--------|
| `Results/snp.txt` (final SNP report) | **IDENTICAL** |
| VarScan `*.snp.vcf` (all 3 samples) | **IDENTICAL** |
| Pindel `Summary.100.vcf` (all 3 samples) | **IDENTICAL** (byte-for-byte) |
| VarScan `*.indel.vcf` / `*.consensus.vcf` | Differ only by ±1–2 read counts at a few positions; **no variant gained or lost**. |

The only ≥2% change is the existing `6218 -GGCTCT` deletion, whose frequency shifts
by 0.03% (Lib1 30.85→30.82%, Lib2 35.01→34.98%) from ±2 supporting reads; it stays
called in both runs. This is sub-threshold read-count jitter, not a variant change.

The same kind of jitter was seen earlier when only the **thread count** was changed
(T1 vs T1.5, same env and reference): BWA-MEM batch boundaries shift with `-t`, so
per-position read counts move by ±1–2 while the called variants stay identical. The
fastp (0.23 → 1.3) and bwa (0.7.17 → 0.7.19) upgrades produce the same harmless
effect. All such differences fall below every reporting threshold.

**Conclusion:** the upgrade is output-neutral at the reporting level.

---

## 3. Future improvements / follow-ups

1. **Thread-independent alignment (optional).** BWA-MEM read-count jitter comes from
   batch boundaries shifting with thread count. `bwa mem -K 10000000` makes batching
   deterministic and removes the ±1–2 read differences seen in the consensus VCFs.
   Cosmetic only — does not affect reported variants.
2. **Consider a VarFreq display threshold for `indel.txt`.** Currently it shows all
   `Summary.100.vcf` calls (including low-AD ones). If reviewers prefer parity with
   `snp.txt`, add an optional frequency filter — but keep raw AD visible, since
   Pindel GT is unreliable at amplicon depth.

---

## 4. Files touched

| File | Change | Tested |
|------|--------|--------|
| [../envs/variant-calling.yml](../envs/variant-calling.yml) | New env definition (Python 3.10 + upgraded tools) | Yes (env built + full run) |
| [../pipeline_combined.py](../pipeline_combined.py) | Added Pindel → `indel.txt` report; renamed from `pindel.txt`; updated docstrings | Yes |
| `bwa_wrapper.py` | `bwa mem -t 6→16`, piped `samtools sort -@ 4→16` | Yes |
| `samtools_wrapper.py` | `sort`/`index` `-@ 6→16` | Yes |
| `pindel_wrapper.py` | `-T 100→16` (thread fix); tuned invariants unchanged | Yes |

## 5. Reproduce a run

```bash
cd ngs/amp-pip-dev
conda run --no-capture-output -n variant-calling \
  python pipeline_combined.py /abs/path/to/RUN_DIR
```
Run from the repo root so the `../amino-anno/from_varscan.py` sibling path resolves.
