# TLA-seq Pipeline — Technical Reference

## 1. Reference construction (`build_reference.py`)

### Host genome

Downloaded once from NCBI using `ncbi-datasets-cli` (`datasets download genome accession`).
The download is unpacked and the first `*_genomic.fna` file found is renamed to `host.fa`.
A GTF is retained if present (optional; used for annotation only). `samtools faidx` is run
to produce the `.fai` index.

**Cache key:** `<host_dir>/<sanitized_accession>/host.fa`  
**Default accession:** `GCF_003668045.3` (Chinese hamster CriGri-PICRH-1.0)  
**Download is atomic:** built into `<acc_dir>.building/` and renamed only on success; an
interrupted download is detected on the next run and retried.

### Mixed reference

Concatenation of `host.fa` and the user-supplied vector FASTA, indexed with `bwa index`
and `samtools faidx`.

**Cache key:** `<mixed_dir>/<accession>__<vector_stem>/mixed.fa`  
**Index completeness check:** `.amb`, `.ann`, `.bwt`, `.pac`, `.sa`, `.fai` must all be
non-empty for a cache hit. Any failure wipes the cache directory so the next run rebuilds.

### Reference path resolution

`pipeline.py` always resolves `host_dir` and `mixed_dir` through `_resolve_reference_dirs`
before calling `build_reference.py`:

| Mode | Trigger | `host_dir` | `mixed_dir` |
|------|---------|-----------|------------|
| Auto (default) | `-r` omitted | `<output>/reference/<accession>/` | `<output>/reference/mixed_reference/` |
| Local | `-r <dir>` | `<dir>/host_genome/` | `<dir>/mixed_ref/` |

The module-level `HOST_GENOME_DIR` / `MIXED_REF_DIR` globals and the `TLA_HOST_GENOME_DIR`
/ `TLA_MIXED_REF_DIR` environment variables are only active when `build_reference.py` is
called standalone.

---

## 2. Preprocessing (`preprocess.py`)

1. **fastp** (v1.3.3): adapter trimming, quality filtering, per-read statistics.
   HTML QC report written to `<sample>/qc/<sample>.fastp.html`.
2. **BFC** (r181): k-mer-based sequencing error correction on the trimmed reads.
   Corrected paired FASTQs feed BWA-MEM.

---

## 3. Alignment (`align.py`)

**BWA-MEM** (0.7.19) against the mixed (host + vector) reference. Supplementary
alignments are kept — split-read segments at integration junctions map to different
loci, so their separate alignments are required for GRIDSS breakend assembly.

Post-alignment: `samtools sort`, `samtools index`, `samtools flagstat` (written to
`<sample>/alignment/<sample>.flagstats.tsv`).

---

## 4. Host coverage and discovery regions

### Coverage track (`host_coverage.py`)

`bedtools genomecov -bg` on the sorted BAM with filter flag `0x104`
(exclude secondary `0x100` and unmapped `0x004`; keep supplementary `0x800` and
duplicate-marked `0x400`). Vector contigs are removed from the bedGraph
(plasmid names read from `plasmid_contigs.txt`).

### Discovery regions (`discovery_sites.py`)

The host bedGraph is binned at **10 Mb** resolution. Active bins are those
containing at least one interval with depth **> 5**. Adjacent active bins and
bins separated by at most **1** inactive bin are merged into candidate regions.
Each region is annotated with `max_depth`, `covered_bp_above_cutoff`, and counts
of GRIDSS breakends, resolved sites, and unresolved junctions that overlap it.

Output: `<sample>/integration/<sample>.discovery_sites.tsv`

---

## 5. GRIDSS2 breakend calling (`gridss_wrapper.py`)

### Deduplication

Before GRIDSS, SAM duplicate-flagged reads (`0x400`) are removed by
`samtools view -F 0x400` into a `.dedup.bam`. This prevents inflated evidence
counts in the host↔host background.

### GRIDSS run

GRIDSS 2.13.2 is run on the dedup BAM with the mixed reference.
Defaults: **8 threads**, **32 g JVM heap**. GRIDSS builds a reference cache
(`mixed.fa.gridsscache`, `mixed.fa.img`) on first use and reuses it.
Both outputs are kept:

- `<sample>_gridss.vcf.gz` — breakend calls (BGZF, tabix-indexed)
- `<sample>_assembly.bam` — GRIDSS-assembled contigs

Skip guard: both files must pass integrity checks (`gzip_ok`, `bam_ok`) for the
GRIDSS step to be skipped; otherwise GRIDSS reruns from its working directory.

---

## 6. Breakend classification (`gridss_sites.py`)

### Single-pass classification

All PASS BND records are read from the compressed VCF as plain text (gzip-readable
BGZF). For each record the mate contig is parsed from the ALT bracket string using
the pattern `[\[\]]([^\[\]]+)\1` to extract `<contig>:<pos>`. Classification is by
whether each anchor contig is in the vector set (from `plasmid_contigs.txt`):

| Class | Condition | Output file |
|-------|-----------|-------------|
| `integration_sites` | one host, one vector anchor | `<sample>.integration_sites.tsv` |
| `vector_sv` | both vector | `<sample>.vector_sv.tsv` |
| `host_host_bnd` | both host | `<sample>.host_host_bnd.tsv` |
| `single_breakend` | no mate bracket in ALT | `<sample>.single_breakend.tsv` |

### EVENT deduplication

Two-sided junctions appear twice in the VCF (a breakend + its `MATEID` partner).
Records are grouped by the `EVENT` INFO tag; each junction is emitted once.

### Host-side anchoring

For `integration_sites`, the row is anchored on the **host side**:
- Prefer the VCF record whose `CHROM` is the host contig (its `CIPOS` is the
  host-breakpoint interval).
- Fallback: if both records have a vector `CHROM`, the host anchor is inferred
  from the mate fields.

### NlaIII annotation (`restriction_site_annotation.py`)

TLA-seq libraries are produced by NlaIII digestion (motif `CATG`, cuts 3' of the
motif on both strands; the motif is self-complementary so one forward scan covers
both strands). After `CATG` is found at position `i` (0-based), the cut position
is stored as `i + 4` (1-based, matching VCF POS convention).

A breakend is "near" a cut site if the distance is ≤ **2 bp**. Each junction is
annotated with:

| Category | Condition |
|----------|-----------|
| `both_ends_near_nlaiii` | both host and vector anchor within 2 bp of a cut |
| `one_end_near_nlaiii` | one anchor within 2 bp |
| `neither_end_near_nlaiii` | neither anchor within 2 bp |

The cut map is stored as `array('I')` for memory efficiency (~5–7× less than a
Python list at whole-genome scale).

---

## 7. Site clustering and resolution (`gridss_split_report.py`)

### Filtering

Junctions with `nlaiii_category == both_ends_near_nlaiii` are excluded from
user-facing reports before clustering. These are likely religation artifacts
(NlaIII self-ligation or ligation of nearby fragments). They remain in the
upstream `integration_sites.tsv`.

### Single-linkage clustering

Reportable junctions are sorted by `(host_contig, host_pos)` and grouped by
single-linkage: junctions within **1 000 bp** of the previous junction on the
same contig are merged into one candidate site. Site IDs are assigned as
`site_001`, `site_002`, …, sorted by `(contig, min_pos, site_id)`.

### Arm orientation

The `arm` column (`left` or `right`) is read from each junction's `breakend`
field: if the ALT string starts with a bracket character (`[` or `]`), the
host-facing arm is `right`; otherwise it is `left`. This preserves the GRIDSS
orientation convention without re-encoding.

### Resolution criterion

A site cluster is **resolved** if it contains at least one `left` arm and at
least one `right` arm. This corresponds to evidence from both sides of the
integration junction. All other reportable junctions are **unresolved**.

Output files:

| File | Content |
|------|---------|
| `<sample>.resolved_sites.tsv` | Junctions from resolved clusters, with `site_id` and `arm` prepended |
| `<sample>.unresolved_junctions.tsv` | Junctions from unresolved clusters |

The two outputs partition all reportable junctions exactly; a `RuntimeError` is
raised if they do not.

---

## 8. Coordinate conventions

- All positions are **1-based** (VCF POS convention), passed through unmodified
  from GRIDSS.
- `cipos` is the GRIDSS `CIPOS` INFO field (`-x,+y`), representing the
  microhomology window around the host breakpoint. Its width equals `homlen`.
- Discovery regions use **0-based half-open** coordinates (bedGraph / bedtools
  convention). Do not mix with VCF positions without adding 1.
- NlaIII cut positions are stored 1-based (offset + motif length = offset + 4).

---

## 9. Vector variant calling (`vector_variants.py`)

VarScan `pileup2snp` and `pileup2indel` are run on an mpileup restricted to the
vector contig (`samtools mpileup -d 0 -A -r <contig>`).

Defaults:

| Parameter | Value |
|-----------|-------|
| Min coverage | 1 500 |
| Min variant allele frequency | 0.01 (1%) |

**Indel type convention:** VarScan encodes deletions as `-BASES` and insertions
as `+BASES` in the `VarAllele` column. The `type` column in the report is
`DEL` when the allele starts with `-` and `INS` when it starts with `+`.
Do not infer type from string lengths.

The `n_pindel_calls` column in all vector characterization reports is always `.`
(Pindel is not wired in the current pipeline).

---

## 10. Two-tier reporting

| Tier | File | Audience |
|------|------|---------|
| Complete evidence | `<sample>.integration_sites.tsv` | Informatics; all GRIDSS PASS host↔vector rows including NlaIII-filtered ones |
| User-facing resolved | `<sample>.resolved_sites.tsv` | Biologist; only paired-arm, non-religation sites |
| User-facing unresolved | `<sample>.unresolved_junctions.tsv` | Biologist; requires IGV review before discarding |

The `final_report/resolved_integration_sites_summary.tsv` concatenates per-sample
resolved sites with a `# Sample: <id>` banner block header. `final_report/integration_summary.tsv`
is a flat table of per-sample counts (one row per sample).

---

## 11. Known limitations and edge cases

- **Multi-copy inserts:** multiple resolved sites at different loci are expected
  and are each reported; no attempt is made to link them as a single event.
- **Complex rearrangements:** inversion or tandem duplication of the vector
  produces vector↔vector GRIDSS breakends in `vector_sv.tsv` but no Phase 2
  realignment or structure inference is performed in v1.
- **Unplaced scaffolds (`NW_` contigs):** treated as host throughout. Coverage
  on unplaced scaffolds can be high without reflecting a true integration; the
  `signal_contigs.png` diagnostic plot helps identify these.
- **Short microhomology:** the CIPOS window can be 0 bp for PRECISE calls at
  blunt junctions. Do not require `homlen > 0` as a quality gate.
- **GRIDSS JVM heap:** 32 g is required for whole-genome CHO assembly. Reducing
  the heap causes silent GRIDSS failures with no VCF output. Increase if the log
  shows Java OutOfMemoryError.
- **`host_host_bnd.tsv`:** kept for manual inspection of TLA proximity-ligation
  junctions and CHO host structural-variant background. Not fed into integration
  calling.
- **Phase 2 not implemented:** per-site realignment of GRIDSS assembly contigs
  and host–vector–host structure strings are not produced in v1.
