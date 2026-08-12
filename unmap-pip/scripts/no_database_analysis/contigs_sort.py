import os
import argparse
from Bio import SeqIO

# Setup command-line argument parsing
parser = argparse.ArgumentParser(description='Sort a FASTA file by sequence length.')
parser.add_argument('fasta', type=str, help='Input FASTA file')
parser.add_argument('-n', type=int, default=None, help='Number of top contigs to include in output')

args = parser.parse_args()

# Read the sequences from the fasta file
sequences = list(SeqIO.parse(args.fasta, "fasta"))

# Sort the sequences by length
sequences.sort(key=lambda seq: -len(seq))

# If -n argument is provided, take only the top n sequencesa
if args.n is not None:
    sequences = sequences[:args.n]

# Get the directory of the input file
input_dir = os.path.dirname(args.fasta)

# Create output file name by appending '.sorted' (and optionally the number of sequences)
output_file_name = os.path.basename(args.fasta).rsplit('.', 1)[0]
if args.n is not None:
    output_file_name += f'.sorted.{args.n}.fasta'
else:
    output_file_name += '.sorted.fasta'

# If the input file has a directory, join it with the output file name
if input_dir:
    output_file = os.path.join(input_dir, output_file_name)
else:
    output_file = output_file_name

# Write the sorted sequences to a new fasta file
SeqIO.write(sequences, output_file, "fasta")
