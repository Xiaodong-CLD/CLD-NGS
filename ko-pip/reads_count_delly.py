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
            svtype = record.info["SVTYPE"]  # Get the SVTYPE (DEL or INS)
            start_pos = record.pos
            end_pos = record.stop
            rr = record.samples[0]['RR']
            rv = record.samples[0]['RV']
            counts.append((start_pos+1, end_pos, rv, svtype,rr))  # Include SVTYPE
    return counts

def write_counts_to_file(output_file_path, total_mapped_reads, delly_counts):
    with open(output_file_path, 'w') as f:
        f.write(f"total number of analyzed reads: {total_mapped_reads}\n")
        f.write("start_pos\tend_pos\treads_count\tSVTYPE\tWT_reads_count\n")  # Add SVTYPE in header
        total_ko_reads = 0
        for count in delly_counts:
            f.write(f"{count[0]}\t{count[1]}\t{count[2]}\t{count[3]}\t{count[4]}\n")  # Write SVTYPE
            total_ko_reads += count[2]
        f.write(f'number of reads without ko: {total_mapped_reads - total_ko_reads}\n')
            


def count_reads_delly(bam_file_path,delly_bcf_file_path,output_file_path):
    total_mapped_reads = get_mapped_reads_count(bam_file_path)
    delly_counts = get_delly_counts(delly_bcf_file_path)
    write_counts_to_file(output_file_path, total_mapped_reads, delly_counts)

if __name__ == "__main__":
    pass

