import argparse
from bwa_wrapper import bwa_runner
from samtools_wrapper import samtools_runner
from varscan_wrapper import varscan_runner
from pindel_wrapper import pindel_runner
from fastp_wrapper import fastp_runner
from plot_coverage_plots import plot_coverage_plots_manager


def quality_check(input_fastq_file_1_path, input_fastq_file_2_path):
    '''
    We are using fastp for adaptive trimming and quality checking of raw fastq files.
    '''
    fastq_file_1_path, fastq_file_2_path = fastp_runner(input_fastq_file_1_path, input_fastq_file_2_path)
    return fastq_file_1_path, fastq_file_2_path


def align_reads(reference_file_path, fastq_file_1_path, fastq_file_2_path, disable_trim_fastq_files):
    '''
    We are using BWA alignment for aligning reads which could be coming from experiments or simulated data.
    We might also use Bowtie in the future.
    '''
    bam_file_path = bwa_runner(reference_file_path, fastq_file_1_path, fastq_file_2_path, disable_trim_fastq_files)
    bam_file_path, sorted_bam_file_path, pileup_file_path = samtools_runner(bam_file_path, reference_file_path)
    return pileup_file_path, sorted_bam_file_path


def call_snvs(pileup_file_path):
    '''
    For SNVs, we are primarily using VarScan. We will also incorporate
    '''
    consensus_file_path=varscan_runner(pileup_file_path)
    return consensus_file_path


def call_indels(sorted_bam_file_path, reference_file_path):
    '''
    We'll be using PINDEL for getting the possible InDels from the BAM files.
    '''
    pindel_runner(sorted_bam_file_path, reference_file_path)

def pipeline_manager():
    #####ARGUMENTS######
    parser = argparse.ArgumentParser()

    parser.add_argument("-r", "--reference-file", help="Path to the reference file.", required=True)
    parser.add_argument("-f1", "--input-file-1", help="Path to input fastq file 1.", required=True)
    parser.add_argument("-f2", "--input-file-2", help="Path to input fastq file 2.", required=True)
    parser.add_argument("-nt", "--disable-trim-fastq-files", help="Disable trimming of fastq files", required=False, action = "store_true")

    args = vars(parser.parse_args())
    reference_file_path = args['reference_file']
    input_fastq_file_1_path = args['input_file_1']
    input_fastq_file_2_path = args['input_file_2']
    disable_trim_fastq_files = args['disable_trim_fastq_files']

    if disable_trim_fastq_files:
        fastq_file_1_path = input_fastq_file_1_path
        fastq_file_2_path = input_fastq_file_2_path
    else:
        #Quality checing of reads and adapter trimming.
        fastq_file_1_path, fastq_file_2_path = quality_check(input_fastq_file_1_path, input_fastq_file_2_path)


    #Align fastq files.
    pileup_file_path, sorted_bam_file_path = align_reads(reference_file_path, fastq_file_1_path, fastq_file_2_path, disable_trim_fastq_files)

    #Call Single Nucleotide Variants.
    consensus_file_path = call_snvs(pileup_file_path)

    #Call InDels.
    call_indels(sorted_bam_file_path, reference_file_path)

    #Create Coverage plot.
    plot_coverage_plots_manager(consensus_file_path, reference_file_path)



if __name__ == "__main__":
    pipeline_manager()
