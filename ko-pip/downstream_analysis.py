"""Per-sample downstream analysis: variant calling + indel/softclip outputs.

Flow (top-to-bottom):

    GATK    ──────► tmp/<sample>/indel/gatk_variants.vcf.gz       (indel caller)
                    tmp/<sample>/indel/gatk_variant.tsv           (per-variant read counts)
    VarScan ──────► tmp/<sample>/snv/                             (SNP + consensus VCFs)
                 └► tmp/<sample>/indel/                           (micro-indel VCF)

    Uninformative composition report (only when GATK called a variant):
        isolate_uninformative_reads   -> Results/<sample>/uninformative_analysis/igv/
        compute_uninformative_summary -> Results/<sample>/uninformative_analysis/

    Breakpoint analysis (always, on the full BAM):
        softclip_analysis           -> Results/<sample>/breakpoint_analysis/
                                       (top_clusters.tsv + top_softclip_sequence.fa;
                                        diagnostic BLAST hits under tmp/<sample>/)

    Unified indel merging (GATK + VarScan2):
        merge_indel_callers         -> Results/<sample>/indel/merged_indel.tsv

    Plots:
        plot_coverage_plots_manager -> Results/<sample>/plots/coverage_plot.jpeg

    The VarScan SNP report (varscan_snp_report) is built from `pipeline.py`
    and written to Results/<sample>/snv/varscan_snp_report.tsv.
"""

import glob
import os
import pathlib
import subprocess

import pysam

from gatk_wrapper import run_gatk, get_mapped_reads_count
from plot_coverage_plots import plot_coverage_plots_manager
from varscan_wrapper import varscan_runner

from merge_indel_callers import merge_indel_callers
from isolate_uninformative_reads import main as isolate_uninformative_reads
from compute_uninformative_summary import main as compute_uninformative_summary
from softclip_analysis import main as softclip_analysis


def _collect_blast_subjects(blast_ref_dir):
    """Return sorted plasmid FASTA paths under blast_ref_dir (empty if None/missing)."""
    if not blast_ref_dir or not os.path.isdir(blast_ref_dir):
        return []
    subjects = []
    for ext in ('*.fa', '*.fasta', '*.fna'):
        subjects += glob.glob(os.path.join(blast_ref_dir, ext))
    return sorted(subjects)


def downstream_analysis(
    reference_path,
    sorted_bam_file_path,
    mpileup_file_path,
    output_dir,
    tmp_sample_dir,
    sample_id=None,
    blast_ref_dir=None,
):
    indel_folder = os.path.join(output_dir, 'indel')
    pathlib.Path(indel_folder).mkdir(parents=True, exist_ok=True)
    # GATK VCFs, gatk_variant.tsv and the VarScan indel VCF are diagnostic
    # intermediates: they live under tmp/<sample>/indel/ so only the merged
    # deliverable (merged_indel.tsv) stays in Results/<sample>/indel/.
    indel_tmp_folder = os.path.join(tmp_sample_dir, 'indel')
    pathlib.Path(indel_tmp_folder).mkdir(parents=True, exist_ok=True)

    # Plasmid FASTAs to BLAST soft clips against (construct-match detection).
    blast_subjects = _collect_blast_subjects(blast_ref_dir)

    # Single canonical sample id reused by every downstream stage so the
    # rows in every output agree.
    _sid = sample_id or pathlib.Path(sorted_bam_file_path).stem.replace('.sorted', '')

    # The whole amplicon is the target: use the first (only) reference contig.
    _chrom = _reference_chrom(reference_path)

    # ------------------------------------------------------------------
    # GATK -> tmp/<sample>/indel/gatk_variants.vcf.gz
    # ------------------------------------------------------------------
    gatk_vcf_path = os.path.join(indel_tmp_folder, 'gatk_variants.vcf.gz')
    if not os.path.exists(gatk_vcf_path):
        try:
            run_gatk(reference_path, sorted_bam_file_path, indel_tmp_folder)
            print('[Takeda Pipeline Message] GATK variant calling executed successfully.')
        except subprocess.CalledProcessError:
            print('[Takeda Pipeline Message] GATK variant calling failed.')

    gatk_found = os.path.exists(gatk_vcf_path)
    if not gatk_found:
        print('[Takeda Pipeline Message] GATK did not find variants; try analyze manually.')

    # ------------------------------------------------------------------
    # Uninformative composition report -- ONLY when GATK called a variant.
    # Isolates the reads GATK ignored at its called sites, summarizes their
    # composition, and writes a BAM + index for IGV. Descriptive only; it
    # does not feed any verdict. Skipped entirely on a GATK no-call (there
    # are no called sites to define "ignored" reads), so the
    # uninformative_analysis/ folder exists iff GATK found a variant.
    # ------------------------------------------------------------------
    if gatk_found:
        try:
            uninf_folder = os.path.join(output_dir, 'uninformative_analysis')
            igv_folder = os.path.join(uninf_folder, 'igv')
            pathlib.Path(igv_folder).mkdir(parents=True, exist_ok=True)

            isolate_uninformative_reads([
                '-b', str(sorted_bam_file_path),
                '-v', str(gatk_vcf_path),
                '-r', str(reference_path),
                '-o', str(igv_folder),
                '--sample', _sid,
            ])

            uninf_bam = os.path.join(igv_folder, f'{_sid}.sorted.bam')

            compute_uninformative_summary([
                '-b', str(uninf_bam),
                '-o', str(os.path.join(uninf_folder, 'uninformative_summary.tsv')),
                '--sample', _sid,
                '--long-softclip-bp', '20',
            ])
            print(f'[Takeda Pipeline Message] Uninformative composition report: {uninf_folder}')
        except Exception as e:
            print(f'[Takeda Pipeline Message] Uninformative composition report failed: {e}')

    # ------------------------------------------------------------------
    # Breakpoint analysis -- ALWAYS run on the full sorted BAM, independent
    # of GATK. Clusters soft-clip junctions and BLASTs the cluster
    # representatives against the plasmid references (construct-match /
    # insert-like flag). Large insertions and construct integrations leave
    # soft-clip junctions here even when no caller emits an indel.
    # ------------------------------------------------------------------
    try:
        bp_folder = os.path.join(output_dir, 'breakpoint_analysis')
        bp_scratch = os.path.join(tmp_sample_dir, 'breakpoint_analysis')
        blast_db_dir = os.path.join(bp_scratch, 'blast_db')
        pathlib.Path(bp_folder).mkdir(parents=True, exist_ok=True)
        pathlib.Path(bp_scratch).mkdir(parents=True, exist_ok=True)

        softclip_args = [
            '-b', str(sorted_bam_file_path),
            '-o', str(bp_folder),
            '--scratch-dir', str(bp_scratch),
            '--blast-db-dir', str(blast_db_dir),
            '--blast-subject', ','.join(blast_subjects),
            '--ref-fasta', str(reference_path),
        ]
        softclip_analysis(softclip_args)

        print(f'[Takeda Pipeline Message] Breakpoint analysis completed: {bp_folder}')
    except Exception as e:
        print(f'[Takeda Pipeline Message] Breakpoint analysis failed: {e}')

    # ------------------------------------------------------------------
    # Call VarScan SNP + indel + consensus VCFs for this sample. The SNP and
    # consensus VCFs are diagnostic intermediates written to tmp/<sample>/snv/;
    # the indel VCF is likewise written to tmp/<sample>/indel/.
    # ------------------------------------------------------------------
    snv_tmp_folder = os.path.join(tmp_sample_dir, 'snv')
    pathlib.Path(snv_tmp_folder).mkdir(parents=True, exist_ok=True)
    output_consensus_vcf_file_path = varscan_runner(
        mpileup_file_path,
        snv_dir=snv_tmp_folder,
        indel_dir=indel_tmp_folder,
    )

    # ------------------------------------------------------------------
    # Merge GATK + VarScan2 indels -> Results/<sample>/indel/merged_indel.tsv
    # ------------------------------------------------------------------
    _varscan_matches = glob.glob(os.path.join(indel_tmp_folder, '*.indel.vcf'))
    varscan_indel_for_merger = _varscan_matches[0] if _varscan_matches else None

    try:
        merged_path = merge_indel_callers(
            gatk_vcf=gatk_vcf_path if gatk_found else None,
            varscan_indel_vcf=varscan_indel_for_merger,
            reference_path=reference_path,
            chrom=_chrom,
            sample_id=_sid,
            output_dir=indel_folder,
            total_mapped_reads=get_mapped_reads_count(sorted_bam_file_path),
        )
        print(f'[Takeda Pipeline Message] Merged indels: {merged_path}')
    except Exception as e:
        print(f'[Takeda Pipeline Message] Indel merger failed: {e}')

    # ------------------------------------------------------------------
    # Plots -> Results/<sample>/plots/
    # ------------------------------------------------------------------
    plots_dir = os.path.join(output_dir, 'plots')
    pathlib.Path(plots_dir).mkdir(parents=True, exist_ok=True)

    # Coverage plot from the VarScan consensus VCF (padded to full amplicon length).
    plot_coverage_plots_manager(output_consensus_vcf_file_path, reference_file_path=reference_path,
                                output_dir=plots_dir)





def _reference_chrom(reference_path):
    """Return the first (and, for an amplicon, only) contig name in the reference."""
    fa = pysam.FastaFile(reference_path)
    try:
        return fa.references[0]
    finally:
        fa.close()


if __name__ == '__main__':
    pass
