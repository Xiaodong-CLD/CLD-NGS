import pysam

def identify_ko_delly(delly_bcf_path, knockout_txt_path):
    with pysam.VariantFile(delly_bcf_path) as bcf, open(knockout_txt_path, "a") as knockout_file:
        for record in bcf:
            # Extract necessary information from the record
            chrom = record.chrom
            pos = record.pos
            end = record.stop
            svtype = record.info["SVTYPE"]
            PE = record.info.get("PE", "NA")
            SR = record.info.get("SR", "NA")
            CIPOS = record.info.get("CIPOS", "NA")
            CIEND = record.info.get("CIEND", "NA")
            
            # Format the information with SVTYPE
            info_string = f"{chrom}\t{pos+1}\t{end}\t{svtype}\n"
            if svtype in ['DEL', 'INS']:
                # Append to the knockout.txt file
                knockout_file.write(info_string)

if __name__ == "__main__":
    pass
