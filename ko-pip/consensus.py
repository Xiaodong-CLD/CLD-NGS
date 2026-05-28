from pyfaidx import Fasta

def generate_consensus(reference_path, ko_path, consensus_path):
    # Load reference sequence
    reference = Fasta(reference_path)
    ref_name = next(iter(reference)).name
    ref_seq = str(next(iter(reference)))

    # Read KO positions from ko_pos.txt
    ko_regions = []
    with open(ko_path, "r") as ko_file:
        for line in ko_file:
            parts = line.strip().split('\t')
            if parts[0] == ref_name:  # ensure the reference name matches
                start, end = int(parts[1]), int(parts[2])
                ko_regions.append((start, end))

    consensus_seqs = [ref_seq]  # Start with the original sequence as a base

    # Generate consensus sequences based on KO regions
    for start, end in ko_regions:
        new_consensus = list(ref_seq)
        new_consensus[start-1:end] = ['-'] * (end-start+1)
        consensus_seqs.append(''.join(new_consensus))

    # Write to the output FASTA file
    with open(consensus_path, "w") as f:
        print(f">original_sequence", file=f)
        for i in range(0, len(ref_seq), 80):
            print(ref_seq[i:i+80], file=f)
        
        for idx, consensus in enumerate(consensus_seqs[1:], 1):
            print(f"\n>consensus_{idx}", file=f)
            for i in range(0, len(consensus), 80):
                print(consensus[i:i+80], file=f)

    print("[Takeda Pipeline Message] Consensus FASTA file has been generated.")

if __name__ == "__main__":
    pass
