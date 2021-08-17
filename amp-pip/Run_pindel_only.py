import argparse
from pindel_wrapper import pindel_runner

def call_indels(sorted_bam_file_path, reference_file_path):
    '''
    We'll be using PINDEL for getting the possible InDels from the BAM files.
    '''
    pindel_runner(sorted_bam_file_path, reference_file_path)

def pipeline_manager():
    #####ARGUMENTS######
    parser = argparse.ArgumentParser()

    parser.add_argument("-r", "--reference-file", help="Path to the reference file.", required=True)
    parser.add_argument("-b", "--sorted-bam-file", help="Path to sorted bam file.", required=True)
    # parser.add_argument("-f2", "--input-file-2", help="Path to input fastq file 2.", required=True)
    # parser.add_argument("-nt", "--disable-trim-fastq-files", help="Disable trimming of fastq files", required=False, action = "store_true")

    args = vars(parser.parse_args())
    reference_file_path = args['reference_file']
    sorted_bam_file_path = args['sorted_bam_file']
    # input_fastq_file_2_path = args['input_file_2']
    # disable_trim_fastq_files = args['disable_trim_fastq_files']

    # if disable_trim_fastq_files:
    #     fastq_file_1_path = input_fastq_file_1_path
    #     fastq_file_2_path = input_fastq_file_2_path
    # else:
    #     #Quality checing of reads and adapter trimming.
    #     fastq_file_1_path, fastq_file_2_path = quality_check(input_fastq_file_1_path, input_fastq_file_2_path)


    #Align fastq files.
    # pileup_file_path, sorted_bam_file_path = align_reads(reference_file_path, fastq_file_1_path, fastq_file_2_path, disable_trim_fastq_files)
    #
    # #Call Single Nucleotide Variants.
    # call_snvs(pileup_file_path)

    #Call InDels.
    call_indels(sorted_bam_file_path, reference_file_path)

    #Create Coverage plot.
    # input_dir = "/".join(input_fastq_file_1_path.split('/')[:-1])
    # plot_coverage_plots_manager(input_dir)



if __name__ == "__main__":
    pipeline_manager()
