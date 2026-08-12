#!/usr/bin/env python3

"""
Pipeline for Knockout analysis (KO-pip) of amplicon sequencing data.

This Python script performs genome assembly and downstream analysis given
a directory of paired-end .fastq files. For each pair of fastq files in
the input directory, it creates a directory in the output folder, where
all the output files from the analysis of that pair will be placed.

The user specifies a directory of paired-end .fastq files (-i), a sample sheet
(.xlsx) as a direct file path (-s), an optional output directory (-o, default:
current directory), a reference FASTA (-r), and a directory of plasmid FASTAs
to BLAST breakpoint clusters against (-br). If the output directory is not
provided, the current working directory is used.

Original script created by Lilin Huang.

Updated by Cindy Ren on 2026-06-22:
- Renamed from pipeline_combined.py; imports updated for the renamed
  preprocess.py / align.py modules.
- All large intermediates (fastp trimmed FASTQ, SAM/BAM, sorted BAM,
  mpileup, depth, reference + index files) are now written under
  <output_dir>/tmp/, with persistent per-sample results kept under
  <output_dir>/Results/<sample>/. Source FASTQ files are read directly
  from the input dir (no more local copy); fastp does the rename to
  <naming>.R1/R2.trimmed.* during trimming.
- fastp HTML/JSON reports are routed to Results/<sample>/qc/ so QC
  deliverables survive even when tmp/ is cleaned.
- Reference FASTA (and any sibling index/.dict files) is copied into
  tmp/reference/ before indexing, so the source reference directory can
  remain read-only.

Per-sample diagnostic intermediates are kept under tmp/<sample>/<module>/ while
Results/<sample>/ holds only the deliverables:
- indel/       : only merged_indel.tsv (GATK/VarScan VCFs live in tmp/)
- consensus/   : only <sample>.consensus.fa (vsearch intermediates + diagnostic
                 tables live in tmp/)
- breakpoint_analysis/ : top_clusters.tsv + top_softclip_sequence.fa (full BLAST
                 hits live in tmp/)
- snv/         : derived deliverables only (varscan_snp_report.tsv; coverage
                 plot under plots/). Raw VarScan SNP/consensus VCFs live in
                 tmp/<sample>/snv/. Amino-acid annotation is archived; the
                 report's AA columns are left blank.
"""
import argparse
import glob
import os
import pandas as pd
import pathlib
import shutil
import subprocess
from preprocess import data_preprocessing
from downstream_analysis import downstream_analysis
from align import assembly
from varscan_snp_report import build_snp_report
from aggregate_run_report import aggregate_run_report
from vsearch_wrapper import vsearch_runner
from gatk_wrapper import get_mapped_reads_count

def name(string):
    '''
    Input: the name of the fatsq files
    Output: the name used to name output files and folders
    '''
    prefix = string.split("_L")[0]
    return prefix

# Parse command line arguments
parser = argparse.ArgumentParser(description='A pipeline for genome assembly.')
parser.add_argument('-i', '--input_dir', required=True, help='Directory containing fastq files')
parser.add_argument('-s', '--sample_sheet', required=True, help='Path to the sample sheet (.xlsx)')
parser.add_argument('-o', '--output_dir', default=os.getcwd(), help='Output directory (default is current directory)')
parser.add_argument('-r', '--reference', required=True, help='Path to the reference genome')
parser.add_argument('-br', '--blast_ref', required=True,
                    help='Directory of plasmid FASTA files (*.fa/*.fasta) to BLAST breakpoint soft-clip clusters against (e.g. data/blast_ref). Required: clustering only groups junctions; BLAST classifies them (construct match / insert like).')
args = parser.parse_args()

# blast_ref is required and must hold at least one plasmid FASTA -- breakpoint
# classification (construct match / insert like) depends on it.
if not os.path.isdir(args.blast_ref):
    parser.error(f"--blast_ref directory not found: {args.blast_ref}")
if not any(glob.glob(os.path.join(args.blast_ref, ext))
           for ext in ('*.fa', '*.fasta', '*.fna')):
    parser.error(f"--blast_ref contains no FASTA (*.fa/*.fasta/*.fna): {args.blast_ref}")

# Validate sample sheet
if not os.path.isfile(args.sample_sheet):
    parser.error(f"--sample_sheet not found: {args.sample_sheet}")
config_f = args.sample_sheet
reference_path = args.reference
sample_records = []


# Read in spreadsheet and convert to dict
config_f_basename = config_f.rstrip('.xlsx')
sheet = pd.read_excel(config_f, index_col=0)
sheet = sheet.fillna('')
d = sheet.to_dict(orient='index')

# Set up a single tmp area to hold all large intermediate files for this
# run: tmp/reference/ for the indexed reference copy, tmp/<sample>/ for
# per-sample FASTQ copies, fastp outputs, SAM/BAM, mpileup, depth, etc.
# Persistent per-sample results stay under Results/<sample>/.
tmp_root = os.path.join(args.output_dir, 'tmp')
tmp_ref_dir = os.path.join(tmp_root, 'reference')
pathlib.Path(tmp_ref_dir).mkdir(parents=True, exist_ok=True)

# Copy the reference FASTA (and any pre-existing sibling index files) into
# tmp/reference/ so we never write to the source reference dir.
ref_basename = os.path.basename(reference_path)
tmp_reference_path = os.path.join(tmp_ref_dir, ref_basename)
if not os.path.exists(tmp_reference_path):
    shutil.copyfile(reference_path, tmp_reference_path)
for _ext in ('.amb', '.ann', '.bwt', '.pac', '.sa', '.fai'):
    _src = reference_path + _ext
    _dst = tmp_reference_path + _ext
    if os.path.exists(_src) and not os.path.exists(_dst):
        shutil.copyfile(_src, _dst)
# GATK .dict lives next to the FASTA with .fa/.fasta replaced by .dict
_ref_dict_src = os.path.splitext(reference_path)[0] + '.dict'
_ref_dict_dst = os.path.splitext(tmp_reference_path)[0] + '.dict'
if os.path.exists(_ref_dict_src) and not os.path.exists(_ref_dict_dst):
    shutil.copyfile(_ref_dict_src, _ref_dict_dst)
reference_path = tmp_reference_path

# Index the (copied) reference if needed
file_format = os.path.splitext(reference_path)[1]
index_file_path = reference_path.replace(file_format, f"{file_format}.amb")
if os.path.exists(index_file_path):
    print("[Takeda Pipeline Message] The reference file has already been indexed.")
else:
    print("[Takeda Pipeline Message] Indexing the reference file.")
    subprocess.run(['bwa', 'index', reference_path], check=True)


folders_list = list()
for sample_idx in d:
    naming = name(d[sample_idx]['fastq (forward)'])
    temp_dir = os.path.join(args.output_dir,"Results")
    sample_dir = os.path.join(temp_dir, naming)
    #make the output directories
    pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)
    # Per-sample tmp area for FASTQ copies + alignment/sort intermediates.
    tmp_sample_dir = os.path.join(tmp_root, naming)
    pathlib.Path(tmp_sample_dir).mkdir(parents=True, exist_ok=True)
    folders_list.append(sample_dir)
    
    fq_clone = d[sample_idx]['Clone ID (optional)']

    fq_1_name = d[sample_idx]['fastq (forward)']
    fq_2_name = d[sample_idx]['fastq (reverse)']
    # Preserve the source extension (.fastq[.gz] or .fq[.gz]) so the trimmed
    # outputs keep the same compression as the input.
    if '.fastq' in fq_1_name:
        ext_base = '.fastq'
    else:
        ext_base = '.fq'
    ext_suffix = fq_1_name.split(ext_base)[1]  # '.gz' or ''
    fq_ext = f"{ext_base}{ext_suffix}"

    # Read source FASTQ files directly — fastp_runner renames + writes
    # trimmed outputs under tmp_sample_dir, and HTML/JSON QC reports under
    # Results/<sample>/qc/ where they persist as deliverables.
    fq_1_path = os.path.join(args.input_dir, fq_1_name)
    fq_2_path = os.path.join(args.input_dir, fq_2_name)
    qc_dir = os.path.join(sample_dir, 'qc')

    # Run fastp_runner + bfc to perform data preprocessing
    processed_fastq_file_1_path, processed_fastq_file_2_path = data_preprocessing(
        fq_1_path, fq_2_path,
        output_dir=tmp_sample_dir,
        naming_prefix=naming,
        ext=fq_ext,
        report_dir=qc_dir,
    )
    if processed_fastq_file_1_path is False:
        raise RuntimeError(
            f"[{naming}] Preprocessing (fastp/bfc) failed — see messages above. "
            f"Check that the FASTQ files are readable: {fq_1_path}"
        )
    
    # Perform assembly (alignment intermediates written under tmp_sample_dir)
    sorted_bam_file_path, mpileup_file_path, _depth_file_path, _collated_bam_file_path = assembly(processed_fastq_file_1_path, processed_fastq_file_2_path,reference_path,tmp_sample_dir)
    if sorted_bam_file_path is False:
        raise RuntimeError(
            f"[{naming}] Assembly (alignment/samtools) failed — see messages above."
        )

    # Perform downstream analysis (persistent results under sample_dir)
        
    downstream_analysis(reference_path, sorted_bam_file_path, mpileup_file_path, sample_dir,
                        tmp_sample_dir,
                        blast_ref_dir=args.blast_ref)
   
     
     
     
    # VSEARCH whole-amplicon consensus haplotypes from the trimmed reads.
    # Persistent deliverable under Results/<sample>/consensus/; vsearch
    # intermediates and diagnostic tables go to tmp/<sample>/consensus/.
    consensus_dir = os.path.join(sample_dir, 'consensus')
    consensus_scratch_dir = os.path.join(tmp_sample_dir, 'consensus')
    pathlib.Path(consensus_dir).mkdir(parents=True, exist_ok=True)
    pathlib.Path(consensus_scratch_dir).mkdir(parents=True, exist_ok=True)
    try:
        consensus_fa = vsearch_runner(
            processed_fastq_file_1_path, processed_fastq_file_2_path,
            reference_path, consensus_dir, consensus_scratch_dir, naming,
        )
        if consensus_fa:
            print(f"[Takeda Pipeline Message] Consensus haplotypes wrote: {consensus_fa}")
        else:
            print("[Takeda Pipeline Message] Consensus (vsearch) produced no haplotypes.")
    except Exception as e:
        print(f"[Takeda Pipeline Message] Consensus (vsearch_wrapper) failed: {e}")

    snp_vcf_folder = os.path.join(sample_dir, 'snv')
    snp_vcf_tmp_folder = os.path.join(tmp_sample_dir, 'snv')
    pathlib.Path(snp_vcf_folder).mkdir(parents=True, exist_ok=True)
    pathlib.Path(snp_vcf_tmp_folder).mkdir(parents=True, exist_ok=True)

    # NOTE: amino-acid annotation is archived. The current annotator assumes the
    # reference is in-frame from position 1, which is wrong for a generic
    # amplicon, so the AA columns were unreliable. A correct future version must
    # take the CDS start position (and strand) to translate variants. Until
    # then, the report is built directly from the plain VarScan SNP table and
    # the amino-acid columns are left blank.
    try:
        snp_report_path = build_snp_report(
            varscan_dir=snp_vcf_tmp_folder,
            sample_id=naming,
            output_dir=snp_vcf_folder,
            total_reads=get_mapped_reads_count(sorted_bam_file_path),
        )
        print(f"[Takeda Pipeline Message] SNP report wrote: {snp_report_path}")
    except Exception as e:
        print(f"[Takeda Pipeline Message] SNP report (varscan_snp_report) failed: {e}")

    sample_records.append((naming, sample_idx, fq_clone))

# Per-run aggregated reports under final_report/.
results_dir = os.path.join(args.output_dir, 'Results')
final_report_dir = os.path.join(args.output_dir, 'final_report')
aggregate_run_report(results_dir, final_report_dir, sample_records)

print(f'[Takeda Pipeline Message] Pipeline complete. {len(sample_records)} sample(s) processed. '
      f'Per-sample results: {results_dir}/ | Aggregated report: {final_report_dir}/')