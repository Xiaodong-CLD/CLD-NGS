import subprocess
import os
from samtools_wrapper import samtools_runner

def alignment(processed_fastq_file_1_path, processed_fastq_file_2_path,reference_path,sample_dir):
    sample_name = os.path.basename(processed_fastq_file_1_path).split(".")[0]
    sam_file_path = os.path.join(sample_dir,f'{sample_name}.sam')
    bam_file_path = os.path.join(sample_dir,f'{sample_name}.bam')
    
    # Check if the output file already exists, if so skip the alignment
    if os.path.exists(bam_file_path):
        print("[Takeda Pipeline Message] Alignment has already been executed. Skip the alignment.")
        return bam_file_path
    try: 
        # Run alignment
        bwa_command = ["bwa", "mem", reference_path, processed_fastq_file_1_path, processed_fastq_file_2_path]
        # Open sam output file
        with open(sam_file_path,"w") as output_file:
            # Run BWA
            subprocess.run(bwa_command, stdout=output_file, check = True)
        
        # Convert sam file to bam file needed for later analysis
        samtools_command = ["samtools","view","-S","-b",sam_file_path]
        # Open bam output file
        with open(bam_file_path, "w") as output_file:
            # Run samtools
            subprocess.run(samtools_command, stdout=output_file, check=True)
        return bam_file_path

    except subprocess.CalledProcessError:
        print("[Takeda Pipeline Message] Alignment failed.")
        return False





def assembly(processed_fastq_file_1_path, processed_fastq_file_2_path,reference_path,sample_dir):
    
    bam_file_path = alignment(processed_fastq_file_1_path, processed_fastq_file_2_path,reference_path,sample_dir)
    sorted_bam_file_path, mpileup_file_path,depth_file_path = samtools_runner(bam_file_path, reference_path)
    
    return sorted_bam_file_path, mpileup_file_path,depth_file_path
    
if __name__ == "__main__":
    pass