import pandas as pd
import os
import glob
import subprocess
import shutil
import pathlib
import argparse
from collections import defaultdict

# Argument parser
parser = argparse.ArgumentParser()
parser.add_argument('input_dir', help='Input directory containing folders')
parser.add_argument('star_combined', help='Path to the directory of combined STAR index')
parser.add_argument('star_human', help='Path to the directory of human STAR index')
parser.add_argument('star_goi', help='Path to the directory of goi STAR index')
parser.add_argument('goi_seq', help='Path to the goi seq file')
args = parser.parse_args()

def name(string):
    prefix = string.split("_L")[0]
    return prefix

def add_newline_if_needed(file_path):
    """Adds a newline at the end of file if it doesn't exist"""
    with open(file_path, 'rb+') as file:
        file.seek(-1, os.SEEK_END)
        last_char = file.read().decode()
        if last_char != '\n':
            print(f"File {file_path} does not end with a newline, adding one...")
            file.write(b'\n')

config_f = glob.glob("*.xlsx")[0]
config_f_basename = config_f.rstrip('.xlsx')
sheet = pd.read_excel(config_f, index_col=0)
sheet = sheet.fillna('')
d = sheet.to_dict(orient='index')

goi_fasta = glob.glob(os.path.join(args.input_dir,'goi','*'))[0]

output_dir = os.path.join(args.input_dir,"Result/")

snp_dfs = []
clone_samples = defaultdict(list)
clone_for_every_sample = True
for sample_idx in d:
    naming = name(d[sample_idx]['fastq (forward)'])
    sample_dir = os.path.join(output_dir, naming)
    pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)
    
    fq_clone = d[sample_idx]['Clone ID (optional)']

    if fq_clone == '':
        clone_for_every_sample = False
    else:
        clone_samples[fq_clone].append(sample_idx)
        fq_new_base = f"{sample_idx}_{fq_clone}"

    # Copy fastq files
    fq_1_name = d[sample_idx]['fastq (forward)']
    if '.fastq' in fq_1_name:
        ext_1 = '.fastq'
    else:
        ext_1 = '.fq'
    fq_1_split = fq_1_name.split(ext_1)
    fq_1_new_name = f"{naming}.R1{ext_1}{fq_1_split[1]}"

    fq_2_name = d[sample_idx]['fastq (reverse)']
    if '.fastq' in fq_2_name:
        ext_2 = '.fastq'
    else:
        ext_2 = '.fq'
    fq_2_split = fq_2_name.split(ext_1)
    fq_2_new_name = f"{naming}.R2{ext_2}{fq_2_split[1]}"

    fq_1_path = os.path.join(args.input_dir, 'Fastq', fq_1_name)
    fq_2_path = os.path.join(args.input_dir, 'Fastq', fq_2_name)
    fq_1_new_path = os.path.join(sample_dir, fq_1_new_name)
    fq_2_new_path = os.path.join(sample_dir, fq_2_new_name)
    shutil.copyfile(fq_1_path, fq_1_new_path)
    shutil.copyfile(fq_2_path, fq_2_new_path)

    fastq_sample_dir = os.path.abspath(sample_dir)

    subprocess.run(["bash", "/mnt/efs/home/zds0975/pipeline_test/voof/RNAseq-var/var_calling.sh", 
                "-d", fastq_sample_dir, 
                "-o", sample_dir,
                "-c", args.star_combined, 
                "-r", args.star_human, 
                "-g", args.star_goi,
                "-s",args.goi_seq], 
                check=True)

    snp_vcf_folder = os.path.join(sample_dir, naming,f"{naming}.snp")
    snp_vcf_path = os.path.join(snp_vcf_folder, f"{naming}.goi.mpileup.snp.vcf")
    snp_vcf_aa_path = os.path.join(snp_vcf_folder, f"{naming}.goi.mpileup.snp.amino_annotated.vcf")
    
     # Run amino acid annotation on varscan snp output
    subprocess.run(["python3",
                    "/mnt/efs/home/zds0975/pipeline_test/voof/amino-anno/from_varscan.py",
                    snp_vcf_path,
                    goi_fasta,
                    "-o",
                    snp_vcf_folder])
    
    snp_df = pd.read_csv(snp_vcf_aa_path, sep='\t')
    snp_df.insert(0, 'SampleID', sample_idx)
    snp_df['CloneID'] = fq_clone
    # Keeping snp output dict for each sample
    snp_dfs.append(snp_df)
    
# Putting all the individual snp dicts together
snp_combined_orig = pd.concat(snp_dfs)
# Dropping unnecessary columns
snp_combined_orig = snp_combined_orig.drop(['Strands1', 'Strands2', 'Qual1', 'Qual2', 'MapQual1', 'MapQual2'], axis=1)
# Converting string percentages to float
snp_combined_orig['VarFreq'] = snp_combined_orig['VarFreq'].apply(lambda x: float(x[:-1]))
# Dropping all entries with VarFreq < 2%
snp_combined_orig = snp_combined_orig[snp_combined_orig['VarFreq'] >= 2]
# Converting back to string
snp_combined_orig['VarFreq'] = snp_combined_orig['VarFreq'].apply(lambda x: str(x) + '%')

# Adding the 'Cons' column to the index
snp_combined = snp_combined_orig.set_index(['SampleID', 'Position', 'Cons'])
snp_combined_orig = snp_combined_orig.drop('CloneID', axis=1)

if clone_for_every_sample:
    by_sample_pos_cons = snp_combined.to_dict('index')
    by_clone_by_pos_cons = defaultdict(lambda: defaultdict(list))
    for k,v in by_sample_pos_cons.items():
        by_clone_by_pos_cons[v['CloneID']][k[1:]].append(k[0])  # k[1:] corresponds to ('Position', 'Cons')

    with open(f"{config_f_basename}.snp.txt", 'w') as out_f:
        for clone_id,sample_ids in clone_samples.items():
            clone_output_df = pd.DataFrame()
            sample_l = clone_samples[clone_id]
            if len(sample_l) > 1:
                pos_cons_samples = by_clone_by_pos_cons[clone_id]
                for pos_cons in pos_cons_samples:
                    if len(pos_cons_samples[pos_cons]) > 1:  # only include if a variant (Position + Cons) appears in more than one sample
                        for s in sample_l:
                            row = snp_combined_orig.query("SampleID == @s and Position == @pos_cons[0] and Cons == @pos_cons[1]")
                            clone_output_df = pd.concat([clone_output_df, row]) 
            else: 
                for s in sample_l:
                        row = snp_combined_orig.query("SampleID == @s")
                        clone_output_df = pd.concat([clone_output_df, row]) 

            out_f.write(f"CloneID: {clone_id}\n")
            out_f.write(f"Replicate Sample IDs: {', '.join(str(x) for x in sample_l)}\n")
            if clone_output_df.empty:
                out_f.write('No positions with variants above 2% frequency in all replicates\n\n')
            else:
                out_f.write(f"{clone_output_df.to_string(index=False)}\n\n")
else:
    with open(f"{config_f_basename}.snp.txt", 'w') as out_f:
        for sample_id in d:
            out_f.write(f"SampleID: {sample_id}\n")
            sample_df = snp_combined_orig.loc[snp_combined_orig['SampleID'] == sample_id] 
            if sample_df.empty:
                out_f.write('No positions with variants above 2% frequency\n\n')
            else:
                out_f.write(f"{sample_df.to_string(index=False)}\n\n")