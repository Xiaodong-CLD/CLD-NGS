# ko-pip Technical Notes

Developer-facing reference for pipeline methods, parameters, formulas, and
implementation decisions.

For installation and execution, see [`README.md`](../README.md).  
For result interpretation, see [`user_guide.md`](user_guide.md).

External output coordinates are 1-based. Internal pysam coordinates are
0-based half-open and converted before writing.

## Pipeline Stages

Entry point: `ko_pip.sh` launches `pipeline.py` in the `ko_pip` conda
environment and a detached `tmux` session.

| Stage | Main module/tool | Main output |
|---|---|---|
| Trim and error-correct | fastp, bfc | Corrected FASTQ and QC reports |
| Align and generate pileup | bwa, samtools | BAM, mpileup, and depth files |
| Call variants | GATK, VarScan2 | Per-sample SNP and indel calls |
| Analyze uninformative reads | pysam modules | `uninformative_analysis/` |
| Analyze soft clips | starcode, blastn | `top_clusters.tsv` |
| Merge indels | bcftools, pysam | `merged_indel.tsv` |
| Build coverage plot | plotly, kaleido | Coverage JPEG |
| Build consensus haplotypes | vsearch | Consensus FASTA |
| Aggregate run | pandas | `final_report/` |

The uninformative-read stage runs only when GATK called at least one variant,
because there are otherwise no called loci at which to classify ignored reads.

If bfc fails, preprocessing continues with fastp-trimmed reads and logs a
warning. This prevents an error-correction failure from discarding an otherwise
usable sample.

## Alignment and Pileup

### Read groups

```text
bwa mem -R '@RG\tID:<sample>\tSM:<sample>\tLB:<sample>\tPL:illumina\tPU:unit1'
```

Read groups are embedded during alignment because GATK requires them, and this
avoids a separate Picard read-group and re-sorting step.

### mpileup settings

```text
samtools mpileup -A -B -d 0
```

- `-A`: retain anomalous or orphan read pairs.
- `-B`: disable BAQ.
- `-d 0`: disable the per-position depth cap.

The depth cap is disabled because the default can truncate high-depth amplicon
evidence and distort VarScan frequencies. BAQ is disabled to avoid suppressing
bases near true indels.

## SNP Calling

VarScan2 command:

```text
varscan mpileup2snp <mpileup> --min-var-freq 0.01
```

No explicit `--min-coverage` is set, so the VarScan2 default applies.

The 1% frequency threshold is intentionally permissive so that low-frequency
signals remain available for review rather than being discarded during calling.

`varscan_snp_report.tsv` is a parsed summary without PASS/FILTER labels.

Additional depth gate:

```text
total_informative >= 0.01 × total_mapped_reads
```

This removes low-coverage primer-end variants that can pass the frequency
threshold despite having little absolute support.

The gate is skipped when total mapped reads are not supplied.

```text
var_freq = alt_support / (ref_support + alt_support)
```

VarScan2 `mpileup2cns` is used for coverage plotting and is not a primary SNP
deliverable.

## Indel Calling and Merging

### GATK HaplotypeCaller

```text
gatk HaplotypeCaller \
  --max-reads-per-alignment-start <2 × mapped_reads> \
  --stand-call-conf 20 \
  --active-probability-threshold 0.0002 \
  --min-base-quality-score 5
```

- The read-start limit effectively disables downsampling so deep amplicon
  evidence is retained.
- `--stand-call-conf 20` sets a moderate reporting-confidence floor.
- The low active-region and base-quality thresholds preserve lower-frequency
  evidence around edited loci.

### VarScan2

```text
varscan mpileup2indel <mpileup> --output-vcf 1
```

`--output-vcf 1` produces a standard VCF that can be normalized with bcftools.
VarScan2 mainly contributes micro-indels and is not relied on for large events.

### Normalization and merge

Both caller VCFs are normalized with:

```text
bcftools norm -f <reference> -m-any
```

Normalization left-aligns alleles and splits multiallelic records so equivalent
calls from the two callers can be merged reliably.

Records merge on normalized `(chrom, start, ref, alt)` and receive:

- `gatk_only`
- `varscan2_only`
- `both`

Spanning-deletion `*` alleles are removed.

### Output coordinates

VCF anchor bases are excluded from the reported event sequence.

For deletions:

```text
start_pos = POS + 1
end_pos   = POS + len(REF) - 1
seq       = REF[1:]
```

For insertions:

```text
start_pos = end_pos = POS
seq       = ALT[1:]
```

This representation is easier to interpret than raw VCF REF/ALT alleles because
it reports only the deleted or inserted sequence.

### Edge and depth gates

```text
EDGE_ARTIFACT_MARGIN_BP = 3
```

`edge_artifact = yes` when an indel starts within 3 bp of either amplicon end.
These calls are commonly caused by primer boundaries or unstable edge alignment.

A `varscan2_only` row is removed when:

```text
AD + RD < 0.01 × total_mapped_reads
```

This removes low-depth VarScan-only calls in samples with much higher overall
coverage and no matching GATK evidence.

Missing caller-side values are written as `.`.

## Soft-Clip Analysis

Reads with flags `4|256|1024|2048` are excluded. These correspond to unmapped,
secondary, duplicate, and supplementary reads.

Leading and trailing soft clips of at least 10 bp are extracted.

```text
CLUSTER_MIN_LEN = 10
```

Shorter clips are usually too nonspecific to support reliable breakpoint
interpretation.

### Starcode clustering

- Algorithm: `--sphere`
- Edit distance: `8`
- Fallback edit distance: `12`
- Top clusters retained: `20`
- Threads: `4`

Clustering is performed separately for `5p` and `3p` clips.

The relatively permissive edit distance collapses sequencing-error variants of
the same junction, while side-specific clustering prevents unrelated left and
right breakpoints from merging.

### BLAST classification

Representative clip sequences are searched against:

- supplied construct/plasmid FASTAs;
- the amplicon reference.

A full hit requires:

```text
percent identity >= 90%
query coverage >= 80%
```

These thresholds require most of the clipped sequence to match at high identity,
rather than classifying short local similarities as full evidence.

Per-sample classifications:

- `construct_blast`: `construct_match`, `partial`, or `.`
- `reference_blast`: `reference_match`, `partial`, or `.`
- per-sample `blast_result`: `hit`, `partial`, or `no hit`

Aggregated `all_breakpoints.tsv` uses:

- `construct`
- `reference`
- `construct,reference`
- `.`

The per-sample and aggregated `blast_result` fields have different value sets
and should not be interpreted interchangeably.

### Large-event flags

`large_insert` requires:

```text
clip length >= 30
support_frac >= 0.01
reference_blast == '.'
```

A long, well-supported non-reference clip is treated as a candidate insertion
or integration even when the inserted sequence cannot be fully resolved.

`large_deletion` requires:

```text
support_frac >= 0.01
reference_blast != '.'
BLAST endpoint at least 5 bp from the anchor
```

The 5 bp distance prevents near-anchor local matches from being labeled as
large deletions.

Both flags may coexist.

Event coordinates are stored in:

```text
DEL;<start>,<end>;<length>
INS;<anchor>,<anchor>;<insert_length>
```

Detailed BLAST coordinates remain in
`tmp/<sample>/breakpoint_analysis/blast_hits_full.tsv` because they are useful
for debugging but too detailed for the main deliverable.

## Consensus Haplotypes

vsearch workflow:

1. Merge paired reads:
   - `--fastq_minovlen 20`
   - `--fastq_maxdiffs 5`
   - `--fastq_maxee 1.0`
   - `--fastq_maxns 0`
2. Dereplicate exact merged reads.
3. Sort by abundance with `--minsize 20`.
4. Retain the top 20 sequences.
5. Align and group sequences into edge-masked haplotypes.

The merge settings require a meaningful overlap and remove merged reads with
excessive disagreement or expected error.

Reporting gates:

```text
consensus_min_size = 50
consensus_min_freq = 0.01
consensus_edge_bp = 5
```

The read-count and frequency gates suppress rare sequence errors, while the edge
mask prevents primer-end instability from defining haplotypes.

Consensus frequency is:

```text
read_support / total_merged_reads
```

Alignment uses affine costs:

- match: `0`
- mismatch: `1`
- gap open: `2`
- gap extension: `1`

The gap-open penalty favors one contiguous indel over several short gaps, which
produces more stable calls in low-complexity sequence.

## Fraction Definitions

All fractions are calculated per sample before run-level filtering. Aggregation
only filters rows, sums breakpoint support, or adds corroboration flags.

| Fraction | Definition |
|---|---|
| `var_freq` | `alt_support / (ref_support + alt_support)` |
| `gatk_alt_frequency_direct` | alternate AD / total AD |
| `gatk_sum_ad_frac` | `sum(AD) / DP_INFO` |
| `varscan2_alt_frequency` | VarScan2 `FREQ` converted to 0–1 |
| `varscan2_sum_ad_frac` | `(AD + RD) / SDP` |
| `support_frac` | cluster reads / analyzed reads |
| `sample_read_frac` | uninformative reads / all sample reads |
| `uninf_alt_frac` | `n_alt / (n_alt + n_ref)` |
| consensus `freq` | haplotype reads / total merged reads |

`gatk_alt_frequency_direct` is preferred over INFO/AF because it is derived
directly from normalized record-level read counts and is easier to audit.

Zero denominators produce `0.0` or `.` rather than an exception.

## Uninformative-Read Accounting

Each GATK-ignored read at a called site is assigned to:

- `ALT`
- `REF`
- `uncountable`

`REF` reads are further divided into `REF_clipped` and `REF_clean`.

`uncountable` is kept as one category because soft-clipped, truncated, and
different-event reads cannot be classified reliably without adding speculative
logic.

Sites remain separate, with one row per GATK alternate allele.

## Run-Level Aggregation

| Output | Keep rule |
|---|---|
| `all_snv.tsv` | `var_freq > 0.01` |
| `all_indel.tsv` | Remove edge artifacts; retain when either caller's attributed-depth fraction is at least `0.01` |
| `all_uninformative.tsv` | `sample_read_frac >= 0.01` |
| `all_breakpoints.tsv` | Merge nearby clusters; retain when summed support is above `0.01` and construct, reference, or large-insert evidence is present |
| `all_consensus.fa` | Retain all reported haplotypes; add `unconfirmed` when no gated indel or breakpoint corroborates the event |

### Breakpoint merging

Breakpoint groups require:

- same side;
- anchors within `3 bp`;
- clip lengths within `3 bp`.

Nearby clusters are merged because sequencing and alignment can shift the same
junction by a few bases.

The representative row contributes `flag` and `variant`; member
`support_frac` values are summed.

### Consensus corroboration

A consensus indel is confirmed when:

- a gated indel of the same type and length lies within `5 bp`; or
- a gated breakpoint lies within `5 bp`.

The 5 bp tolerance absorbs left-alignment and anchor shifts in low-complexity
regions.

Only DEL and INS portions of `sig` are checked. SNV and MNV portions are not
used for corroboration.

## Key Constants

| Parameter | Value |
|---|---|
| SNP minimum variant frequency | `0.01` |
| SNP informative-depth fraction | `0.01` |
| VarScan-only indel depth fraction | `0.01` |
| Run-level attributed-depth gates | `0.01` |
| Breakpoint support gates | `0.01` |
| Uninformative-read gate | `0.01` |
| Consensus minimum frequency | `0.01` |
| GATK stand-call confidence | `20` |
| GATK active probability | `0.0002` |
| GATK minimum base quality | `5` |
| Edge-artifact margin | `3 bp` |
| Breakpoint position tolerance | `3 bp` |
| Breakpoint length tolerance | `3 bp` |
| Consensus confirmation tolerance | `5 bp` |
| Large-deletion anchor tolerance | `5 bp` |
| Minimum soft-clip length | `10 bp` |
| Large-insert clip length | `30 bp` |
| Starcode edit distance | `8`, fallback `12` |
| Consensus minimum reads | `50` |

## Assumptions and Limitations

- Single-contig amplicon reference.
- No amino-acid annotation without explicit CDS start and strand.
- Uninformative-read analysis requires a GATK call.
- Large events rely primarily on soft-clip and consensus evidence.
- Breakpoint coordinates are estimates derived from anchors and BLAST endpoints.
- fastp currently uses unpinned default trimming behavior.
- Automated gates support prioritization but do not replace review of
  borderline calls.

## Items to Revisit

- Consolidate repeated `0.01` thresholds only if they are intended to change 
  together.
- Review whether the separate `3 bp` tolerances should remain independent.
- Document fastp settings if custom trimming thresholds are introduced.