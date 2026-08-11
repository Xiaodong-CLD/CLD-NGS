# ko-pip: Knockout Amplicon Pipeline

`ko-pip` analyzes paired-end amplicon sequencing from knockout experiments. It reports:

- **SNVs** with VarScan2
- **Small and medium indels** with GATK and VarScan2
- **Large deletions, insertions, and construct integrations** from soft-clipped reads
- **Whole-amplicon consensus haplotypes**

Results are produced per sample and summarized across the full run.

For result interpretation, see [docs/user_guide.md](docs/user_guide.md).  
For methods, parameters, and formulas, see [docs/technical_notes.md](docs/technical_notes.md).

## Requirements

- Linux
- `conda`
- `tmux`
- Conda environment created from [`envs/ko_pip.yml`](envs/ko_pip.yml)

Key dependencies include Python 3.10, fastp, bfc, bwa, samtools, bcftools,
GATK4, VarScan2, vsearch, starcode, BLAST, pysam, pandas, plotly, and kaleido.

## Inputs

Required inputs:

| Input | Description |
|---|---|
| **FASTQ directory** (`-i`) | Directory containing the paired-end `.fastq[.gz]` / `.fq[.gz]` files |
| **Sample sheet** (`-p`) | Path to the `.xlsx` sample sheet (direct file path, not a directory) |
| **Reference FASTA** (`-r`) | Single-contig amplicon reference |
| **BLAST reference directory** (`-b`) | Directory of construct/plasmid FASTAs (`*.fa`, `*.fasta`, or `*.fna`) |

Required sample-sheet columns:

| Column | Description |
|---|---|
| `fastq (forward)` | R1 filename inside the FASTQ directory |
| `fastq (reverse)` | R2 filename inside the FASTQ directory |
| `Clone ID (optional)` | Optional clone label |

FASTQs may be compressed or uncompressed: `.fastq`, `.fastq.gz`, `.fq`, or `.fq.gz`.

## Run

```bash
./ko_pip.sh \
  -i <fastq_dir> \
  -p <sample_sheet.xlsx> \
  -o <output_dir> \
  -r <reference.fa> \
  -b <blast_ref_dir> \
  [-s <tmux_session>]
```

| Flag | Description |
|---|---|
| `-i` | Directory containing paired FASTQ files |
| `-p` | Path to the sample sheet (`.xlsx`) |
| `-o` | Output directory |
| `-r` | Amplicon reference FASTA |
| `-b` | Construct/plasmid BLAST-reference directory |
| `-s` | Optional tmux session name |

The pipeline runs in a detached `tmux` session.

```bash
tail -f logs/<session>.log
tmux attach -t <session>
```

## Output Structure

```text
<output_dir>/
  tmp/                         # intermediate and diagnostic files
    reference/                 # indexed reference copy
    <sample>/                  # per-sample intermediates (BAM, VCF, etc.)
  Results/
    <sample>/
      qc/
      snv/
        varscan_snp_report.tsv
      indel/
        merged_indel.tsv
      consensus/
        <sample>.consensus.fa
      breakpoint_analysis/
        top_clusters.tsv
        top_softclip_sequence.fa
      uninformative_analysis/  # only when GATK called a variant
        uninformative_summary.tsv
        igv/
      plots/
        coverage_plot.jpeg
  final_report/
    all_snv.tsv
    all_indel.tsv
    all_breakpoints.tsv
    all_uninformative.tsv
    all_consensus.fa
    coverage_plots/
```

Files under `Results/` and `final_report/` are deliverables. Files under `tmp/`
are retained for debugging and may be deleted after successful review.

## Main Deliverables

**Run-level** (`final_report/`):

| File | Purpose |
|---|---|
| `final_report/all_indel.tsv` | Primary small/medium indel summary |
| `final_report/all_breakpoints.tsv` | Large deletion, insertion, and construct-junction evidence |
| `final_report/all_consensus.fa` | Whole-amplicon haplotypes |
| `final_report/all_snv.tsv` | SNVs above the reporting threshold |
| `final_report/all_uninformative.tsv` | Reads ignored by GATK at called sites |
| `final_report/coverage_plots/` | Per-sample coverage plots |

**Per-sample** (`Results/<sample>/`):

| File | Purpose |
|---|---|
| `indel/merged_indel.tsv` | Merged GATK + VarScan indel calls for this sample |
| `snv/varscan_snp_report.tsv` | VarScan SNV calls for this sample |
| `consensus/<sample>.consensus.fa` | Whole-amplicon consensus haplotypes |
| `breakpoint_analysis/top_clusters.tsv` | Top soft-clip breakpoint clusters |
| `breakpoint_analysis/top_softclip_sequence.fa` | Soft-clip cluster sequences (BLAST input) |
| `uninformative_analysis/uninformative_summary.tsv` | Uninformative-read summary at GATK-called sites |
| `uninformative_analysis/igv/` | BAM files of uninformative reads for IGV inspection |
| `plots/coverage_plot.jpeg` | Per-sample coverage plot |

Run-level TSVs use `# sample:` banner blocks to separate samples.

## Recommended Review Order

1. `all_indel.tsv`
2. `all_breakpoints.tsv`
3. `all_consensus.fa`
4. `all_snv.tsv`
5. `all_uninformative.tsv`, when caller support is unexpectedly low

See the [user guide](docs/user_guide.md) for interpretation examples and
confidence guidance.

## Limitations

- Supports a single-contig amplicon reference only.
- Does not report amino-acid consequences.
- Large events may require breakpoint evidence and orthogonal validation.
- The uninformative-read report is generated only when GATK called a variant.
- Automated filters support review but do not replace inspection of borderline cases.
