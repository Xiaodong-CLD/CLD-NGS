#!/usr/bin/env python

'''
This script is a modification of the script located at:
    https://github.com/ShivamSharma13/genome-assembly-pipeline/blob/master/spades_wrapper.py

The above mentioned script on github is written by me (Shivam Sharma).
'''
import subprocess
import os
import glob
from filter_deletions import filter_deletions_from_pindel2vcf


def get_pindel_intermittent_files(pindel_output_file_prefix_path):
    #The following suffixes for Long Indels, Short Indels, Deletions and Inversions are the only files we are considering in the pipeline.
    suffixes_of_interest = ['_SI', '_LI', '_D']
    return [pindel_output_file_prefix_path + i for i in suffixes_of_interest if os.path.exists(pindel_output_file_prefix_path + i)]


def pindel_runner(sorted_bam_file_path, reference_file_path):
    #Check if PINDEL direcotry exists, if not, create one.
    if not os.path.exists(os.path.join("/".join(sorted_bam_file_path.split('/')[:-1]), 'pindel-output')):
        os.mkdir(os.path.join("/".join(sorted_bam_file_path.split('/')[:-1]), 'pindel-output'))

    config_file_path = os.path.join("/".join(sorted_bam_file_path.split('/')[:-1]), 'pindel-output', sorted_bam_file_path.split('/')[-1].replace(".bam", ".pindel.config.txt"))
    pindel_output_file_prefix_path =  os.path.join("/".join(sorted_bam_file_path.split('/')[:-1]), 'pindel-output', sorted_bam_file_path.split('/')[-1].replace(".bam", ".pindel"))

    #Since PINDEL is very sensitive, we will use several parameters to get more confident results.
    #Always make sure the string starts and ends with a space.
    pindel_filter_query = " -E 0.9 "

    #Create a config file.
    config_file_text = sorted_bam_file_path + '\t' + '398' + '\t' + sorted_bam_file_path.split('/')[-1]
    with open(config_file_path, 'w') as f:
        f.write(config_file_text)

    #Check if PINDEL intermittent files have already been produced.
    #[NOTE] PINDEL intermittent files are created when PINDEL starts running but the tool takes a lot of time to run.
    #If unsure if the tool ran successfully, please delete the intermittent files produced by PINDEL and run the pipeline again.
    if not os.path.exists(pindel_output_file_prefix_path + '_SI'):
        #Running PINDEL to produce the intermittent files.
        try:
            print("[Takeda Pipeline Message] Running PINDEL on sample {}".format(sorted_bam_file_path.split('/')[-1]))
            print("[Takeda Pipeline Message] PINDEL will produce intermittent files for various structural variants such as: SI, LI, INV, D etc...")
            bash_output = subprocess.call("pindel -f " + reference_file_path + " -i " + config_file_path + pindel_filter_query + " -T 16 -o " + pindel_output_file_prefix_path, shell = True)

        except subprocess.CalledProcessError:
            print("[Takeda Pipeline Message] PINDEL failed running. :( ")
            return False
    else:
        print("[Takeda Pipeline Message] Found PINDEL intermittent files, using them...")


    #Running pindel2vcf to produce the vcf files from the detailed intermittent files produced by PINDEL.
    input_file_dir_path = "/".join(sorted_bam_file_path.split('/')[:-1])
    pindel_output_vcf_files_dir_path = os.path.join(input_file_dir_path, 'pindel-output')

    #Check output directory for PINDEL VCF files.
    if not os.path.exists(pindel_output_vcf_files_dir_path):
        print("[Takeda Pipeline Message] Creating subdirectory for PINDEL vcf files...")
        os.mkdir(pindel_output_vcf_files_dir_path)

    ###############################################################################################################
    ############################################Creating PINDEL VCF files.#########################################
    ###############################################################################################################

    pindel_intermittent_files_paths = get_pindel_intermittent_files(pindel_output_file_prefix_path)

    # Run pindel2vcf once at -e 1 (ceiling); fractional AD filtering is applied
    # downstream in pipeline_combined.parse_pindel_summary (ad_frac >= 1%).
    for pindel_intermittent_file_path in pindel_intermittent_files_paths:
        print("[Takeda Pipeline Message] Generating vcf files from PINDEL intermittent files...")
        pindel_output_vcf_file_path = os.path.join(
            pindel_output_vcf_files_dir_path,
            pindel_intermittent_file_path.split('/')[-1] + '.1.vcf')
        print(pindel_output_vcf_file_path)

        #Check if pindel2vcf files already exist.
        if os.path.exists(pindel_output_vcf_file_path):
            print("[Takeda Pipeline Message] VCF file for PINDEL output {} already exist, skipping...".format(pindel_output_vcf_file_path.split('/')[-1]))
            continue
        try:
            bash_output = subprocess.call("pindel2vcf -p " + pindel_intermittent_file_path + " -r " + reference_file_path + " -R x -d 00000000  -e 1 -v " + pindel_output_vcf_file_path, shell = True)
        except subprocess.CalledProcessError:
            print("[Takeda Pipeline Message] PINDEL2VCF failed... :( ")
            return False

        #Deletion files generated by PINDEL have a lot of RPL type InDels and those which have HOMLEN & HOMSEQ.
        #Right now we will filter those out and only carry forqard with the remaining Deletions.
        if "_D" in pindel_intermittent_file_path:
            filter_deletions_from_pindel2vcf(pindel_output_vcf_file_path)

    sample_name = pindel_output_vcf_files_dir_path.split('/')[-1]
    _build_summary_vcf(pindel_output_vcf_files_dir_path, sample_name)


def _build_summary_vcf(vcf_dir, sample_name):
    # vcf_dir (str): the directory where the all the vcf files are saved,
    #                should be something like '/home/adam/' or '/home/adam'
    head_set = []
    SI_infors, D_filter_infors, LI_infors = [],[],[]
    for files_path in glob.glob(os.path.join(vcf_dir,"*.1.vcf")):
        # Walk through the files in vcf_directory to find files that
        # match our needs. Target_files are ended with 'D_filter.1.vcf',
        # or 'SI.1.vcf' or 'LI.1.vcf'
        if files_path.endswith("SI.1.vcf"):
            head_set, SI_infors = obtain_infor_from_file('### SI: Short Insertion', files_path, head_set)
        elif files_path.endswith("D_filter.1.vcf"):
            head_set, D_filter_infors = obtain_infor_from_file('### D_filter: Filtered Delection', files_path, head_set)
        elif files_path.endswith("LI.1.vcf"):
            head_set, LI_infors = obtain_infor_from_file('### LI: Long Insertion', files_path, head_set)

    combined_infors = head_set + SI_infors + LI_infors + D_filter_infors

    with open(os.path.join(vcf_dir,"{}.Summary.vcf".format(sample_name)),'w') as file_handle:
        for item in combined_infors:
            file_handle.write("{}\n".format(item))


def obtain_infor_from_file(starter_infor, file_path, head_set):
    file_related_infors = [starter_infor]
    with open(file_path, "r") as file_handle:
        lines = file_handle.readlines()
        for index, line in enumerate(lines):
            stripped_line = line.strip()
            if stripped_line.startswith('##'):
                if stripped_line not in head_set:
                    ## meta data lines, only append on its first occurence
                    head_set.append(stripped_line)
            else:
                file_related_infors.append(stripped_line)
    return head_set, file_related_infors


if __name__ == "__main__":
    #This script will be called by the main script.
    pass
