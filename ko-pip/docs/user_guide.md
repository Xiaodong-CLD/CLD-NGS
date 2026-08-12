# ko-pip User Guide

This guide explains which outputs to review and how to interpret them. It is
intended for scientific users rather than pipeline developers.

## Start with the Run-Level Report

Open `final_report/` and review files in this order:

1. **`all_indel.tsv`** — primary knockout call
2. **`all_breakpoints.tsv`** — large events missed by standard callers
3. **`all_consensus.fa`** — independent whole-amplicon haplotypes
4. **`all_snv.tsv`** — additional substitutions
5. **`all_uninformative.tsv`** — troubleshooting weak caller evidence

Coverage plots are under `final_report/coverage_plots/`.

## 1. Indel Results

`all_indel.tsv` combines GATK and VarScan2 calls.

Important columns:

| Column | Meaning |
|---|---|
| `start_pos`, `end_pos` | 1-based event coordinates |
| `seq` | Deleted or inserted sequence |
| `type` | `DEL` or `INS` |
| `length` | Event length |
| `source` | `gatk_only`, `varscan2_only`, or `both` |
| `gatk_alt_frequency_direct` | GATK alternate-allele fraction |
| `varscan2_alt_frequency` | VarScan2 alternate-allele fraction |
| `gatk_sum_ad_frac` | Fraction of raw locus depth assigned by GATK |
| `varscan2_sum_ad_frac` | Fraction of raw locus depth assigned by VarScan2 |

Missing caller values are written as `.`.

### Interpreting allele fractions

| Observed fraction | Often consistent with |
|---|---|
| Near 1.0 | Predominant edited allele, homozygous/hemizygous edit, or allele imbalance |
| Near 0.5 | Heterozygous edit or two similarly represented alleles |
| Below 0.2 | Low-level signal, mosaicism, or technical artifact |

Allele fractions should be interpreted together with caller agreement,
attributed-depth fractions, consensus haplotypes, and breakpoint evidence.

### Practical confidence guide

| Evidence | Interpretation |
|---|---|
| Both callers agree, strong fraction, high attributed depth | High confidence |
| One caller only, but supported by consensus | Moderate confidence |
| Low fraction or low attributed depth | Review supporting outputs |
| Event occurs near an amplicon edge | Likely edge/primer artifact |

Events marked `edge_artifact = yes` are excluded from the aggregated
`all_indel.tsv`.

Large deletions and insertions may not appear in this file. Continue to
`all_breakpoints.tsv` even when `all_indel.tsv` is empty.

## 2. Breakpoint Results

`all_breakpoints.tsv` summarizes clusters of soft-clipped reads. These reads
contain sequence that could not be aligned continuously to the amplicon.

Important columns:

| Column | Meaning |
|---|---|
| `side` | `5p` or `3p` side of the aligned read |
| `anchor` | 1-based junction position on the amplicon |
| `summed_support_frac` | Combined fraction of reads supporting the locus |
| `blast_result` | `construct`, `reference`, `construct,reference`, or `.` |
| `flag` | `large_deletion`, `large_insert`, both, or `.` |
| `variant` | Estimated event coordinates |

Interpretation:

- **`blast_result = construct`**: clipped sequence matches a supplied
  construct/plasmid reference.
- **`blast_result = reference`**: clipped sequence matches another region of
  the amplicon and may indicate a large deletion.
- **`flag = large_deletion`**: breakpoint geometry supports a deletion.
- **`flag = large_insert`**: a long, well-supported non-reference clip supports
  an insertion or integration.

For a flagged large deletion, check whether `all_indel.tsv` contains a matching
deletion. Use caller-derived allele frequency when the same event is resolved by
a variant caller.

A construct match or unresolved large event should usually be confirmed with
an orthogonal method such as targeted PCR or Sanger sequencing.

## 3. Consensus Haplotypes

`all_consensus.fa` reports abundant full-length amplicon sequences reconstructed
independently of GATK and VarScan2.

Example header:

```text
>hap1; read_support=194506; freq=0.900691; covered=2-130; n_internal_variants=1; insertions=none; sig=71-78:DEL:CTTTGTTA; flags=none
```

Key fields:

| Field | Meaning |
|---|---|
| `read_support` | Number of merged reads assigned to the haplotype |
| `freq` | Fraction of merged reads assigned to the haplotype |
| `covered` | Reference span represented |
| `sig` | Compact variant description |
| `insertions` | Reported insertions |
| `flags` | Warnings such as `low_support` or `unconfirmed` |

`flags=unconfirmed` means an indel in the haplotype was not corroborated by a
gate-passing indel or breakpoint result. Treat it cautiously.

The current `sig` field reports deletions; small insertions are listed under
`insertions=` and are not cross-checked against callers.

## 4. SNV Results

`all_snv.tsv` contains SNVs with `var_freq > 0.01`.

Important columns:

| Column | Meaning |
|---|---|
| `position`, `ref`, `alt` | SNV position and alleles |
| `ref_support`, `alt_support` | Informative read counts |
| `var_freq` | Alternate-allele fraction |
| `p_value` | VarScan2 p-value |
| `alt_strands_plus`, `alt_strands_minus` | Alternate reads by strand |

The `consensus` column uses IUPAC nucleotide codes. The report does not assign
PASS/FILTER labels or amino-acid consequences.

Low-frequency SNVs should not be interpreted as knockout evidence without
supporting biological context.

## 5. Uninformative-Read Results

`all_uninformative.tsv` summarizes reads that overlap a GATK-called site but
were not used by GATK.

Important columns:

| Column | Meaning |
|---|---|
| `sample_read_frac` | Ignored reads at the site divided by all sample reads |
| `n_alt` | Ignored reads supporting the alternate event |
| `n_ref` | Ignored reads supporting the reference |
| `n_uncountable` | Reads that could not be classified |
| `uninf_alt_frac` | `n_alt / (n_alt + n_ref)` |

Use this file when attributed-depth fractions in `all_indel.tsv` are unexpectedly
low.

- `uninf_alt_frac` near 0: ignored reads mostly support reference.
- Higher `uninf_alt_frac`: GATK may have undercounted alternate-supporting reads.

This report is generated only for samples where GATK called at least one
variant.

## Overall Call Guidance

| Evidence combination | Suggested interpretation |
|---|---|
| Both callers agree and consensus confirms | High-confidence edited allele |
| One caller plus consensus support | Moderate confidence |
| Caller evidence without consensus support | Review coverage and sequence context |
| Breakpoint-only large deletion or insertion | Candidate large event; validate |
| Construct match without a resolved insertion call | Candidate construct integration; validate |
| Only low-frequency SNVs | Negative, low-level, or ambiguous result |

Automated categories are screening aids. Borderline or biologically important
results should be reviewed manually and, when needed, validated experimentally.

## Suggested Screenshots for a Visual Guide

A Word or PDF version can use the same content with:

1. A screenshot of `final_report/`
2. One clear indel example from `all_indel.tsv`
3. The corresponding coverage plot
4. The matching consensus haplotype
5. One breakpoint-only or ambiguous example

Keep the Markdown guide as the version-controlled source and update the visual
copy when a new pipeline version is released.
