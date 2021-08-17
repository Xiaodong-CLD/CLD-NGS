import os
import argparse
import subprocess
from pipeline import pipeline_manager
from merge_fastq import gather_fastq_file_pairs

'''
Pipeline runner is small wrapper that can be used to run the [Takeda Pipeline] on multiple samples in one go.
One sample should consist 2 fastq files [Paired End].

Expected input for this script is a directory which can have multiple dierectories containing paired fastq data.
'''
def run_pipeline_on_fastq_files(input_directory, reference_file_path, disable_trim_fastq_files):
    '''
    We will use os.walk() to get all the subdirectoris and get the files from each of those directories.
    '''
    fastq_files = {}
    for path, directories, _ in os.walk(input_directory):
        for directory in directories:
            fastq_file_pairs = gather_fastq_file_pairs(os.path.join(path, directory))
            for fastq_file_1, fastq_file_2 in fastq_file_pairs.items():
                #Running the pipeline now.
                if "fastq" in fastq_file_1 or "fq" in fastq_file_1:
                    if "sai" not in fastq_file_1 and "trimmed" not in fastq_file_1:
                        fastq_file_1_path = os.path.join(path, directory, fastq_file_1)
                        fastq_file_2_path = os.path.join(path, directory, fastq_file_2)     
          
                        try:
                            print("[Takeda Pipeline Message] Now running variant calling pipeline for fastq files: {}, {}".format(fastq_file_1, fastq_file_2))
                            if disable_trim_fastq_files:
                                bash_output = subprocess.call("python3 pipeline.py -r " + reference_file_path + " -f1 " + fastq_file_1_path + " -f2 " + fastq_file_2_path + " -nt ", shell = True)
                            else:
                                bash_output = subprocess.call("python3 pipeline.py -r " + reference_file_path + " -f1 " + fastq_file_1_path + " -f2 " + fastq_file_2_path, shell = True)                        

                        except subprocess.CalledProcessError:
                            print("Running pipeline failed.")
    return True
            

def pipeline_runner():
    #####ARGUMENTS######   
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--input-directory", help="Path to input directory which will have subdirectories with sample fastq files.", required=True)
    parser.add_argument("-r", "--reference-file", help="Path to reference file.", required=True)
    parser.add_argument("-nt", "--disable-trim-fastq-files", help="Disable trimming of fastq files.", required=False, action='store_true')

    args = vars(parser.parse_args())
    input_directory_path = args['input_directory']
    reference_file_path = args['reference_file']
    disable_trim_fastq_files = args['disable_trim_fastq_files']

    '''Get the fastq files and run the pipeline.'''
    run_pipeline_on_fastq_files(input_directory_path, reference_file_path, disable_trim_fastq_files)
    

if __name__ == "__main__":
    pipeline_runner()
