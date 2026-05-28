import pysam

def get_mapped_reads_count(bam_file_path):
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")
    mapped_reads = bamfile.mapped
    bamfile.close()
    return mapped_reads

def get_gatk_counts(vcf_file_path):
    counts = []
    with pysam.VariantFile(vcf_file_path) as vcf:
        for record in vcf:
            # Determine if it's a deletion or an insertion
            svtype = "DEL" if len(record.ref) > len(record.alts[0]) else "INS"

            start_pos = record.pos
            end_pos = record.pos + len(record.ref) - 1  # Adjusted for insertions
            ad_values = record.samples[0]['AD']
            reads_count = ad_values[1] if len(ad_values) > 1 else 0
            wt_reads_count = ad_values[0] if len(ad_values) > 1 else 0
            counts.append((start_pos+1, end_pos, reads_count, svtype,wt_reads_count))
    return counts


def write_counts_to_file(output_file_path, total_mapped_reads, gatk_counts):
    with open(output_file_path, 'w') as f:
        f.write(f"total number of analyzed reads: {total_mapped_reads}\n")
        f.write("start_pos\tend_pos\treads_count\tSVTYPE\tWT_reads_count\n")  # Add SVTYPE in header
        total_ko_reads = 0
        for count in gatk_counts:
            f.write(f"{count[0]}\t{count[1]}\t{count[2]}\t{count[3]}\t{count[4]}\n")  # Write SVTYPE
            total_ko_reads += count[2]
        f.write(f'number of reads without ko: {total_mapped_reads - total_ko_reads}\n')

def count_reads_gatk(bam_file_path, gatk_vcf_file_path, output_file_path):
    total_mapped_reads = get_mapped_reads_count(bam_file_path)
    gatk_counts = get_gatk_counts(gatk_vcf_file_path)
    write_counts_to_file(output_file_path, total_mapped_reads, gatk_counts)

if __name__ == "__main__":
    pass
