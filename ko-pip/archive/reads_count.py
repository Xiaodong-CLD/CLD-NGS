import pysam


def get_mapped_reads_count(bam_file_path):
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")
    mapped_reads = bamfile.mapped
    bamfile.close()
    return mapped_reads

def get_delly_counts(bcf_file_path):
    counts = []
    with pysam.VariantFile(bcf_file_path) as bcf:
        for record in bcf:
            start_pos = record.pos
            end_pos = record.stop
            rr = record.samples[0]['RR']
            rv = record.samples[0]['RV']
            reads_count = rr + rv
            counts.append((start_pos, end_pos, reads_count))
    return counts

def write_counts_to_file(output_file_path, total_mapped_reads, delly_counts):
    with open(output_file_path, 'w') as f:
        f.write(f"total number of analyzed reads: {total_mapped_reads}\n")
        f.write("start_pos\tend_pos\treads_count\n")
        for count in delly_counts:
            f.write(f"{count[0]}\t{count[1]}\t{count[2]}\n")

def count_reads(bam_file_path,delly_bcf_file_path,output_file_path):
    total_mapped_reads = get_mapped_reads_count(bam_file_path)
    delly_counts = get_delly_counts(delly_bcf_file_path)
    write_counts_to_file(output_file_path, total_mapped_reads, delly_counts)

if __name__ == "__main__":
    pass

