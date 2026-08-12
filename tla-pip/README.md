# tla-pip — TLA-seq Integration-Site Pipeline

Maps transgene integration sites in single-cell clones from paired-end
TLA-seq (Targeted Locus Amplification) FASTQ files.

## Workflow

```
FASTQ pairs
  → fastp trimming + BFC error correction
  → BWA-MEM alignment to mixed (host + vector) reference
  → GRIDSS2 breakend calling
  → NlaIII annotation and artifact filtering
  → host-vector junction clustering and resolution
  → VarScan vector variant calling
  → per-run aggregate report
```

## Requirements

Conda environment `tla_pip` (see `envs/tla-pip.yml`). Key tools:
`bwa 0.7.19`, `samtools 1.23.1`, `fastp 1.3.3`, `bfc r181`,
`bedtools 2.31.1`, `gridss 2.13.2` (requires `openjdk`),
`varscan 2.4.6`, `ncbi-datasets-cli 16.41.0`, `tmux`.

## Inputs

| Flag | Required | Description |
|------|----------|-------------|
| `-s` | yes | Sample sheet (`.xlsx` or `.tsv`) |
| `-f` | yes | Directory containing FASTQ files |
| `-v` | yes | Vector/GOI FASTA (shared across all samples in the run) |
| `-o` | yes | Output directory |
| `-a` | no | Host genome NCBI accession (default: `GCF_003668045.3` — CriGri-PICRH-1.0) |
| `-r` | no | Local reference directory (development override; see below) |
| `--force-download` | no | Re-download host genome even if cached |

**Sample sheet required columns:** `Sample ID (unique)`, `fastq (forward)`, `fastq (reverse)`

## References

The `data/` directory is **not tracked in git** (large files). To reproduce
this run, place the following files in the paths shown before running the
pipeline with `-r data/`.

### Host genome — CHO CriGri-PICRH-1.0

| Field | Value |
|---|---|
| Species | *Cricetulus griseus* (Chinese hamster), strain 17A/GY |
| Assembly | **CriGri-PICRH-1.0** |
| NCBI accession | **GCF_003668045.3** |
| Download | https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003668045.3/ |
| md5 (genomic.fna) | `0ab5dfbed6405f6b09f0c2e9d183f9c5` |
| md5 (genomic.gtf) | `808c19e9f4ef6534e735d77391156579` |

Expected path after download: `data/host_genome/GCF_003668045.3/`

```bash
# Download via NCBI Datasets CLI (included in tla_pip conda env)
datasets download genome accession GCF_003668045.3 \
    --include genome,gtf --filename ncbi_dataset.zip
unzip ncbi_dataset.zip -d data/host_genome/GCF_003668045.3/
```

### Vector / GOI FASTA

| Field | Value |
|---|---|
| Contig name | `D08_HI_HS_AvB1` |
| Length | 9 906 bp |
| Description | Knock-in construct used for TLA-seq experiments in this repo |
| File | Provide the original FASTA as `-v vector.fa` at run time |

The pipeline concatenates this with the host genome to build the mixed
reference. The mixed reference cache is stored at
`data/mixed_ref/<accession>__<vector_contig>/`.

### Reference download (auto mode)

By default the host genome is downloaded from NCBI once into
`<output>/reference/<accession>/` and reused. The mixed reference
(host + vector) is cached under `<output>/reference/mixed_reference/`.

To use a pre-built reference (development, CI):

```bash
-r /path/to/refs   # must contain host_genome/ and mixed_ref/ subdirs
```

## Running

```bash
# Create and activate the environment once
conda env create -f envs/tla-pip.yml

# Production run (reference auto-downloaded into output/)
bash tla_pip.sh -s samples.xlsx -f Fastq/ -v vector.fa -o runs/my_run/

# Development run with local reference
bash tla_pip.sh -s samples.xlsx -f Fastq/ -v vector.fa \
    -o runs/dev_run/ -r /shared/refs/

# Force re-download
bash tla_pip.sh -s samples.xlsx -f Fastq/ -v vector.fa \
    -o runs/my_run/ --force-download
```

Log: `<output>/logs/<session>.log`. Attach: `tmux attach -t <session>`.

## Main outputs

```
<output>/
  logs/                          pipeline log
  reference/                     host genome + mixed reference caches
  final_report/
    integration_summary.tsv      per-sample breakpoint/site counts
    resolved_integration_sites_summary.tsv
    vector_snp_summary.tsv
    vector_indel_summary.tsv
    vector_sv_summary.tsv
    coverage_plots/
  <sample_id>/
    integration/
      <sample>.resolved_sites.tsv
      <sample>.unresolved_junctions.tsv
      <sample>.discovery_sites.tsv
      <sample>.summary.tsv
    vector_characterization/
      <sample>.*.snv_report.tsv
      <sample>.*.indel_report.tsv
      <sample>.vector_structural_breakpoints.tsv
    coverage/
      <sample>.host_coverage.png
      <sample>.vector_coverage.png
    alignment/
      <sample>.flagstats.tsv
    qc/
      <sample>.fastp.html
```

## Documentation

- [User Guide](docs/user_guide.md) — inputs, running, interpreting results, troubleshooting
- [Technical Reference](docs/technical.md) — algorithms, thresholds, coordinate conventions

## Limitations

- Designed for clonal TLA-seq libraries; multi-copy or mosaic inserts may appear as multiple sites.
- Single-vector-contig assumption (one FASTA per run). Samples with different vectors must be run separately.
- Phase 2 (per-site realignment and host–vector–host structure consensus) is not yet implemented.
- The `n_pindel_calls` column in vector characterization reports is always `.` (Pindel not wired).
