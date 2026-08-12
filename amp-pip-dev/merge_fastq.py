#!/mnt/efs/home/cui3008/miniconda3/envs/variant-calling/bin/python
'''
This script takes an input directory and merges all the fastq files (paired) 
present and outputs a final fastq file.
'''
import argparse
import os

def gather_fastq_file_pairs(input_dir):
    '''
    This function goes to a directory and returns all the fastq file pairs present in it.
    Expected output will be something like:
        {'a1.fq': 'a2.fq', 'b1.fq':'b2.fq'}
    '''
    files = [i for i in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, i))]

    #Get first (_1) fastq files.
    first_files = [i for i in files if "_1" in i and  ".sai" not in i]
    unique_file_ids = [i.replace("_1.fastq", "") for i in first_files]

    return {i+"_1.fastq":i+"_2.fastq" for i in unique_file_ids}


def create_output_fastq_files(fastq_file_map, output_file_path, input_dir):
    output_file_1_path = output_file_path + "_1.fastq"
    output_file_2_path = output_file_path + "_2.fastq"

    with open(output_file_1_path, 'ab') as w1, open(output_file_2_path, 'ab') as w2:
        for file_1, file_2 in fastq_file_map.items():
            with open(os.path.join(input_dir, file_1), 'rb') as f:
                w1.write(f.read())
            with open(os.path.join(input_dir, file_2), 'rb') as f:
                w2.write(f.read())

    return True


def merge_reads_manager():
    #Take in the reference file as arguments.
    parser = argparse.ArgumentParser()

    #Add argument options
    parser.add_argument("-i", "--input_dir", help="Path of directory which has fastq files (paired).", required = True)
    parser.add_argument("-o", "--output_file_path", help="Output file path. Output files will be printed with suffixes _1 & _2. For example: '-o output/ten_concentration'", required = True)

    #Fetch the argument data.
    args = vars(parser.parse_args())
    input_dir = args['input_dir']
    output_file_path = args['output_file_path']

    #Get mapped fastq files.
    fastq_file_map = gather_fastq_file_pairs(input_dir)
    create_output_fastq_files(fastq_file_map, output_file_path, input_dir)


if __name__ == "__main__":
    merge_reads_manager()
