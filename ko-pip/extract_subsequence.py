import sys
from Bio import SeqIO

def extract_subsequence(reference_fa, start, end, output_file):
    for record in SeqIO.parse(reference_fa, "fasta"):
        SeqIO.write(record[start:end], output_file, "fasta")

if __name__ == "__main__":
    reference_fa = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    output_file = sys.argv[4]
    extract_subsequence(reference_fa, start, end, output_file)
