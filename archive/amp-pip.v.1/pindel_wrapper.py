#!/usr/bin/env python

'''
This script is a modification of the script located at: 
	https://github.com/ShivamSharma13/genome-assembly-pipeline/blob/master/spades_wrapper.py

The above mentioned script on github is written by me (Shivam Sharma).
'''
import subprocess
import os
import random
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
    #pindel2vcf_filter_query = " -e 30 "

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
            bash_output = subprocess.call("pindel -f " + reference_file_path + " -i " + config_file_path + pindel_filter_query + " -T 100 -o " + pindel_output_file_prefix_path, shell = True)

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

    pindel2vcf_filter_queries = [' -e 15 ', ' -e 30 ', ' -e 50 ', ' -e 75 ', ' -e 100 ']
    pindel_intermittent_files_paths = get_pindel_intermittent_files(pindel_output_file_prefix_path)



    for pindel2vcf_filter_query in pindel2vcf_filter_queries:
        for pindel_intermittent_file_path in pindel_intermittent_files_paths:
            print("[Takeda Pipeline Message] Generating vcf files from PINDEL intermittent files...")
            pindel_output_vcf_file_path = os.path.join(pindel_output_vcf_files_dir_path, pindel_intermittent_file_path.split('/')[-1] + "."  + pindel2vcf_filter_query.replace("-e", "").replace(" ", "") + '.vcf')
            print(pindel_output_vcf_file_path)

            #Check if pindel2vcf files already exist.
            if os.path.exists(pindel_output_vcf_file_path):
                print("[Takeda Pipeline Message] VCF file for PINDEL output {} already exist, skipping...".format(pindel_output_vcf_file_path.split('/')[-1]))
                continue
            try:
                bash_output = subprocess.call("pindel2vcf -p " + pindel_intermittent_file_path + " -r " + reference_file_path + " -R x -d 00000000 " + pindel2vcf_filter_query + "-v " + pindel_output_vcf_file_path, shell = True)
            except subprocess.CalledProcessError:
                print("[Takeda Pipeline Message] PINDEL2VCF failed... :( ")
                return False

            #Deletion files generated by PINDEL have a lot of RPL type InDels and those which have HOMLEN & HOMSEQ.
            #Right now we will filter those out and only carry forqard with the remaining Deletions.
            if "_D" in pindel_intermittent_file_path:
                filter_deletions_from_pindel2vcf(pindel_output_vcf_file_path) 

            
if __name__ == "__main__":
	#This script will be called by the main script.
	pass
