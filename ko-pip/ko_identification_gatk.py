import pysam

def identify_ko_gatk(gatk_vcf_path, knockout_txt_path):
    with pysam.VariantFile(gatk_vcf_path) as vcf, open(knockout_txt_path, "a") as knockout_file:
        for record in vcf:
            chrom = record.chrom
            start_pos = record.pos
            svtype = "NA"  # SVTYPE is not directly available in VCF, so we deduce it

            if len(record.ref) > len(record.alts[0]):  # This indicates a deletion
                end_pos = record.pos + len(record.ref) - 1
                svtype = "DEL"
            elif len(record.ref) < len(record.alts[0]):  # This indicates an insertion
                end_pos = record.pos  # For insertions, the end position is typically the start position
                svtype = "INS"

            if svtype in ["DEL", "INS"]:
                info_string = f"{chrom}\t{start_pos+1}\t{end_pos}\t{svtype}\n"
                knockout_file.write(info_string)

if __name__ == "__main__":
    pass
