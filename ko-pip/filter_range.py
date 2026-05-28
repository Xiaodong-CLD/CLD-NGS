import os

def read_ranges(non_knockout_file):
    ranges = []
    with open(non_knockout_file, 'r') as f:
        for line in f:
            segments = line.strip().split('\t')
            for segment in segments:
                start, end = map(int, segment.split('-'))
                ranges.append((start, end))
    return ranges

def filter_vcf(vcf_file, ranges):
    filtered_records = []
    with open(vcf_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            try:
                position = int(fields[1])
            except ValueError:  # It's likely a header or invalid line
                filtered_records.append(line)
                continue
            
            if any(start <= position <= end for start, end in ranges):
                filtered_records.append(line)
    return filtered_records


def write_to_file(filtered_records, output_file):
    with open(output_file, 'w') as f:
        for record in filtered_records:
            f.write(record)

def filter_based_on_ranges(non_knockout_file, vcf_file):
    ranges = read_ranges(non_knockout_file)
    filtered_records = filter_vcf(vcf_file, ranges)
    output_file = vcf_file.replace(".vcf", "_filtered.vcf")
    
    # Check if the filtered file already existed
    if os.path.exists(output_file):
        print(f"[Takeda Pipeline Message] snp.vcf file has already been created in {output_file}")
        return output_file
    
    write_to_file(filtered_records, output_file)
    return output_file

if __name__ == "__main__":
    pass
