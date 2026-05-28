import subprocess
import os
from varscan_wrapper import varscan_runner
from plot_coverage_plots import plot_coverage_plots_manager
from consensus import generate_consensus
from ko_identification_delly import identify_ko_delly
from ko_identification_gatk import identify_ko_gatk
from non_knockout import extract_non_knockout_intervals
from reads_count_delly import count_reads_delly
from reads_count_gatk import count_reads_gatk
import pathlib

def downstream_analysis(processed_fastq_file_1_path,processed_fastq_file_2_path,reference_path,sorted_bam_file_path, mpileup_file_path,depth_file_path,output_dir):

    current_dir = pathlib.Path(__file__).parent.resolve()
###################################################################################################################################       
    delly_folder = os.path.join(output_dir,'delly')
    # Run variant calling and deletion finding
    vd_folder = delly_folder
    if os.path.exists(vd_folder):
        print(f"[Takeda Pipeline Message] The directory {vd_folder} already exists.")
    else:
        print(f"[Takeda Pipeline Message] Create output directory: {vd_folder}.")
        os.makedirs(vd_folder)
    delly_path = os.path.join(vd_folder,"delly.bcf")
    variants_path = os.path.join(vd_folder,"variants.vcf.gz")
    if os.path.exists(delly_path) and os.path.exists(variants_path):
        print(f"[Takeda Pipeline Message] The variants and deletion calling already performed.")
    else:
        try:
            script_path = os.path.join(current_dir, 'variant_deletion_calling.sh')
            subprocess.check_call([script_path,reference_path,sorted_bam_file_path,vd_folder])
            print("[Takeda Pipeline Message] The variants and deletion calling executed successfully.")
        except subprocess.CalledProcessError:
            print("[Takeda Pipeline Message] The variants and deletion calling failed") 
  
# ******************************************************************************************** 
    delly_found_del = True
    if os.path.exists(delly_path):
        analysis_output = os.path.join(output_dir, 'analysis')
        analysis_output = os.path.join(output_dir, 'delly')
        if os.path.exists(analysis_output):
            print(f"[Takeda Pipeline Message] The directory {analysis_output} already exists.")
        else:
            print(f"[Takeda Pipeline Message] Create output directory: {analysis_output}.")
            os.makedirs(analysis_output)
    
        # Write the knock-out position file
        knockout_file = os.path.join(analysis_output,"knockout_pos.txt")
        if os.path.exists(knockout_file):
            print(f"[Takeda Pipeline Message] The file {knockout_file} already exists. Skipping generating knockout position file.")
        else:
            print(f"[Takeda Pipeline Message] Generating the file for knockout position, output file is: {knockout_file}.")
            identify_ko_delly(delly_path,knockout_file)
        # Write the non-knockout file
        non_knockout_file = os.path.join(analysis_output,"non_knockout_pos.txt")
        if os.path.exists(non_knockout_file):
            print(f"[Takeda Pipeline Message] The file {non_knockout_file} already exists. Skipping generating non_knockout position file.")
        else:
            print(f"[Takeda Pipeline Message] Generating the file for non_knockout position, output file is: {non_knockout_file}.")
            extract_non_knockout_intervals(reference_path, knockout_file, non_knockout_file) 
            
        # Write the reads count file
        reads_count_file = os.path.join(analysis_output,"reads_count.txt")
        if os.path.exists(reads_count_file):
            print(f"[Takeda Pipeline Message] The file {reads_count_file} already exists. Skipping generating reads count file.")
        else:
            print(f"[Takeda Pipeline Message] Generating the reads count file, output file is: {reads_count_file}.")
            count_reads_delly(sorted_bam_file_path, delly_path, reads_count_file)
        
        # Generate consensus reads
        consensus_path = os.path.join(analysis_output,'consensus.fasta')
        
        if os.path.exists(consensus_path):
            print(f"[Takeda Pipeline Message] The file {consensus_path} already exists.")
        else:
            print(f"[Takeda Pipeline Message] Create file: {consensus_path}.")
            generate_consensus(reference_path,knockout_file, consensus_path)
    else:
        print(f"[Takeda Pipeline Message] Delly does not found variants check gatk results.")
        result_txt = os.path.join(delly_folder, 'delly_not_found_ko.txt')
        with open(result_txt,'w') as f:
            f.write("Delly does not found variants.")
        delly_found_del = False
            

# ******************************************************************************************** 
    # GATK Variant Calling
    gatk_folder = os.path.join(output_dir, 'gatk')
    gatk_vcf_path = os.path.join(gatk_folder, "gatk_variants.vcf.gz")
    if not os.path.exists(gatk_vcf_path):
        try:
            script_path_gatk = os.path.join(current_dir, 'gatk_variant_calling.sh')
            subprocess.check_call([script_path_gatk, reference_path, sorted_bam_file_path, gatk_folder])
            print("[Takeda Pipeline Message] GATK variant calling executed successfully.")
        except subprocess.CalledProcessError:
            print("[Takeda Pipeline Message] GATK variant calling failed")
    if os.path.exists(gatk_vcf_path):
        analysis_output = os.path.join(output_dir, 'analysis')
        analysis_output = os.path.join(output_dir, 'gatk')
        if os.path.exists(analysis_output):
            print(f"[Takeda Pipeline Message] The directory {analysis_output} already exists.")
        else:
            print(f"[Takeda Pipeline Message] Create output directory: {analysis_output}.")
            os.makedirs(analysis_output)
    
        # Write the knock-out position file
        knockout_file = os.path.join(analysis_output,"knockout_pos.txt")
        if os.path.exists(knockout_file):
            print(f"[Takeda Pipeline Message] The file {knockout_file} already exists. Skipping generating knockout position file.")
        else:
            print(f"[Takeda Pipeline Message] Generating the file for knockout position, output file is: {knockout_file}.")
            identify_ko_gatk(gatk_vcf_path,knockout_file)
        # Write the non-knockout file
        non_knockout_file = os.path.join(analysis_output,"non_knockout_pos.txt")
        if os.path.exists(non_knockout_file):
            print(f"[Takeda Pipeline Message] The file {non_knockout_file} already exists. Skipping generating non_knockout position file.")
        else:
            print(f"[Takeda Pipeline Message] Generating the file for non_knockout position, output file is: {non_knockout_file}.")
            extract_non_knockout_intervals(reference_path, knockout_file, non_knockout_file) 
            
        # Write the reads count file
        reads_count_file = os.path.join(analysis_output,"reads_count.txt")
        if os.path.exists(reads_count_file):
            print(f"[Takeda Pipeline Message] The file {reads_count_file} already exists. Skipping generating reads count file.")
        else:
            print(f"[Takeda Pipeline Message] Generating the reads count file, output file is: {reads_count_file}.")
            count_reads_gatk(sorted_bam_file_path, gatk_vcf_path, reads_count_file)
        
        # Generate consensus reads
        consensus_path = os.path.join(analysis_output,'consensus.fasta')
        
        if os.path.exists(consensus_path):
            print(f"[Takeda Pipeline Message] The file {consensus_path} already exists.")
        else:
            print(f"[Takeda Pipeline Message] Create file: {consensus_path}.")
            generate_consensus(reference_path, knockout_file, consensus_path)        


    else:
        result_txt = os.path.join(gatk_folder, 'gatk_not_found_ko.txt')
        with open(result_txt,'w') as f:
            f.write("GATK does not found variants.")
        if delly_found_del == False:
            print("[Takeda Pipeline Message] Both Delly and GATK not found deletion, try analyze manually.")
            analysis_output = os.path.join(output_dir, 'analysis')
            result_txt = os.path.join(analysis_output, 'not_fount_ko.txt')
            with open(result_txt,'w') as f:
                f.write("Both Delly and GATK not found deletion, try analyze manually.")
        
    
        
#####################################################################################################################################

   
    # Variant Calling
    output_consensus_vcf_file_path = varscan_runner(mpileup_file_path)
    
    # Generate Plot
    plot_coverage_plots_manager(output_consensus_vcf_file_path)
    
    return non_knockout_file
    


if __name__ == "__main__":
    pass