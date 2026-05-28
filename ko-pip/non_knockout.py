from pyfaidx import Fasta

def extract_non_knockout_intervals(reference_path, knockout_path, output_path):
    # Load the reference fasta
    reference = Fasta(reference_path)
    ref_length = len(str(reference[list(reference.keys())[0]]))

    with open(output_path, 'w') as out_file:
        with open(knockout_path, 'r') as f:
            for line in f:
                _, start, end, *_ = line.strip().split()
                start, end = int(start), int(end)

                non_knockouts_for_allele = []
                
                if start > 1:
                    non_knockouts_for_allele.append(f"{1}-{start-1}")
                
                if end < ref_length:
                    non_knockouts_for_allele.append(f"{end+1}-{ref_length}")
                
                out_file.write("\t".join(non_knockouts_for_allele) + "\n")

if __name__ == "__main__":
    pass
