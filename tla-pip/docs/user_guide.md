# TLA-seq Pipeline — User Guide

## 1. Preparing inputs

### Sample sheet

A `.xlsx` or tab-separated `.tsv` file with exactly these column headers:

| Column | Description |
|--------|-------------|
| `Sample ID (unique)` | Alphanumeric identifier, no spaces |
| `fastq (forward)` | Filename (not full path) of the R1 FASTQ |
| `fastq (reverse)` | Filename of the R2 FASTQ |

All FASTQ files listed must exist in the directory passed to `-f`. Additional
columns are ignored.

### Vector FASTA

A single FASTA file containing the vector/GOI sequence. All samples in a run
share the same file (passed with `-v`). Samples with different vectors must be
run separately.

---

## 2. Production run (automatic reference download)

```bash
conda activate tla_pip
bash tla_pip.sh \
    -s samples.xlsx \
    -f /path/to/Fastq/ \
    -v vector.fa \
    -o /path/to/runs/my_run/
```

On first run the host genome (`GCF_003668045.3` — Chinese hamster CriGri-PICRH-1.0)
is downloaded from NCBI into `<output>/reference/` and cached. Subsequent runs
reuse it. Estimated run time: ~2–4 hours per sample depending on library depth
(GRIDSS dominates).

**Follow progress:**
```bash
tail -f /path/to/runs/my_run/logs/my_run.log
# or attach to the tmux session
tmux attach -t my_run
```

To use a different host genome accession:
```bash
bash tla_pip.sh ... -a GCF_XXXXXXXXX.X
```

To force re-download even if a cache exists:
```bash
bash tla_pip.sh ... --force-download
```

---

## 3. Development / test run (local reference)

If a reference has already been downloaded and indexed in a shared location,
pass it with `-r`. The directory must already exist and must contain
`host_genome/` and `mixed_ref/` subdirectories (the layout that a previous
auto-download produced).

```bash
bash tla_pip.sh \
    -s samples.xlsx \
    -f Fastq/ \
    -v vector.fa \
    -o runs/dev/ \
    -r /shared/refs/
```

No NCBI download is attempted when `-r` is supplied.

---

## 4. Checking a successful run

After completion, verify:

1. Log ends with `TLA pre-analysis pipeline complete.`
2. `final_report/integration_summary.tsv` has one data row per sample with
   non-zero `n_host_vector_breakpoints`.
3. Each `<sample>/integration/<sample>.summary.tsv` exists and has one data row.
4. Coverage plots appear in `final_report/coverage_plots/`.

If a sample is skipped, search the log for `skipping sample` and the error line
above it.

---

## 5. Interpreting results

### 5.1 Integration summary (`final_report/integration_summary.tsv`)

One row per sample. Key columns:

| Column | Meaning |
|--------|---------|
| `n_host_vector_breakpoints` | Raw GRIDSS PASS host↔vector breakends |
| `n_host_vector_both_ends_nlaiii` | Subset with both ends near an NlaIII site — likely religation artifacts (excluded from site reports) |
| `n_host_vector_breakpoints_retained` | Breakpoints carried into site calling after NlaIII filtering |
| `n_resolved_sites` | Integration sites with evidence from both host-flanking arms |
| `n_unresolved_junctions` | Remaining host-vector junctions without a matching paired arm |
| `n_vector_vector_breakpoints` | Within-vector GRIDSS breakends (SV candidates) |

### 5.2 Resolved integration sites (`<sample>/integration/<sample>.resolved_sites.tsv`)

Each row is one junction arm; rows sharing a `site_id` belong to the same
candidate site. A resolved site has both a `left` and a `right` arm.

Key columns:

| Column | Meaning |
|--------|---------|
| `site_id` | Cluster identifier (e.g. `site_001`) |
| `arm` | `left` or `right` (host-facing orientation from GRIDSS ALT bracket) |
| `host_contig` | Host chromosome |
| `host_pos` | Host breakpoint position (1-based, VCF POS) |
| `cipos` | GRIDSS CIPOS interval — the microhomology window |
| `homlen` / `homseq` | Microhomology length and sequence |
| `vector_contig` / `vector_pos` | Vector-side anchor |
| `breakend` | Raw GRIDSS ALT string encoding orientation |
| `qual` / `vf` / `sr` / `rp` | GRIDSS quality, variant fraction, split-read, read-pair support |
| `nlaiii_category` | `one_end_near_nlaiii` or `neither_end_near_nlaiii` |

**Highest-confidence sites:** PRECISE (no `IMPRECISE` flag), high `qual`, `vf > 0`, both arms present, `one_end_near_nlaiii`.

### 5.3 Unresolved junctions (`<sample>/integration/<sample>.unresolved_junctions.tsv`)

Same schema as resolved sites but the cluster has only one arm. These are real
GRIDSS breakends that passed NlaIII filtering but could not be paired — either
the second arm is below GRIDSS detection threshold, was truncated, or is in a
complex rearrangement. Do not discard; review IGV before concluding the site is
absent.

### 5.4 Discovery regions (`<sample>/integration/<sample>.discovery_sites.tsv`)

Coarse host-enrichment regions from binned coverage (10 Mb bins, cutoff depth > 5).
Annotated with the number of GRIDSS breakends, resolved sites, and unresolved
junctions falling within each region. Used as an independent sanity check:
every resolved site should fall within a discovery region. A discovery region
with no GRIDSS evidence but high coverage warrants manual review.

### 5.5 Vector variants

- `<sample>/vector_characterization/<sample>.<contig>.snv_report.tsv` — VarScan SNVs (min coverage 1 500, min VAF 1%)
- `<sample>/vector_characterization/<sample>.<contig>.indel_report.tsv` — VarScan indels; `type` is `DEL` for VarScan `-BASES` alleles and `INS` for `+BASES`
- `<sample>/vector_characterization/<sample>.vector_structural_breakpoints.tsv` — GRIDSS vector↔vector breakends

### 5.6 Coverage plots

- `<sample>/coverage/<sample>.host_coverage.png` — whole-genome host coverage in 100 kb bins (lane plot)
- `<sample>/coverage/<sample>.vector_coverage.png` — per-base vector coverage (linear, unannotated)

The aggregate copies live in `final_report/coverage_plots/` prefixed with the
sample ID.

---

## 6. When to use IGV or wet-lab validation

- Any `resolved` site with `imprecise = True` or `qual < 100` (⚠ threshold not hard-coded; inspect distribution in your run)
- All `unresolved` junctions before concluding a site is absent
- `n_resolved_sites > 2` for a clonal sample (may indicate a mixed population)
- Any vector variant with VAF between 10–90% (suggests partial integration or mixed clone)
- A discovery region with no matched GRIDSS evidence

---

## 7. Common errors

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `FileNotFoundError: GOI/vector FASTA not found` | `-v` path wrong or file empty | Check the path; run `wc -c vector.fa` |
| `ValueError: --references-dir does not exist` | `-r` path typo or not yet created | Create the directory or omit `-r` |
| `No genome FASTA (*.fna) found in NCBI download` | NCBI download incomplete or network error | Delete `<output>/reference/<accession>.building` and rerun |
| `Alignment failed for <sample>` | BWA index missing from mixed ref cache | Delete `<output>/reference/mixed_reference/` and rerun |
| `Integration-site identification failed` | GRIDSS out-of-memory | Increase JVM heap: edit `DEFAULT_JVMHEAP` in `gridss_wrapper.py` (default `32g`) |
| Sample skipped, no `*.sorted.bam` | Preprocessing failed | Check `<output>/<sample>/qc/<sample>.fastp.html` and the log |
| `conda env 'tla_pip' not found` | Environment not yet created | Run `conda env create -f envs/tla-pip.yml` |

---

## 8. Worked interpretation example

Given a sample with:
- `n_host_vector_breakpoints = 4`, `n_host_vector_both_ends_nlaiii = 2`
- `n_host_vector_breakpoints_retained = 2`, `n_resolved_sites = 1`, `n_unresolved_junctions = 0`

Interpretation: 2 raw breakends were filtered as likely religation artifacts (both NlaIII ends). The remaining 2 form one resolved site (`site_001`) with a left and a right arm. This is consistent with a single-copy integration at one locus. Open `<sample>.resolved_sites.tsv`, confirm both rows share `site_id = site_001`, note `host_contig` and `host_pos`, and cross-check that `host_pos` falls within a discovery region in `<sample>.discovery_sites.tsv`. Check `homlen` — a value of 4–5 bp of microhomology at a TTAA/TTTAA motif is typical for this library preparation.
