#!/usr/bin/env python3

import argparse
from collections import defaultdict
import glob
import os
import pandas as pd
import pathlib
import shlex
import shutil
import subprocess

# Wrapper for amplicon pipeline that calls pipeline.py.
# Supports running amplicon pipeline from a configuration spreadsheet,
# and takes replicates into account when generating a report.
# Includes the amino acid annotation on the VarScan snp report.

def name(string):
    prefix = string.split("_L")[0]
    return prefix

current_dir = pathlib.Path(__file__).parent.resolve()

parser = argparse.ArgumentParser()
parser.add_argument('input_dir', help='Path to input directory')
args = parser.parse_args()

snp_dfs = []
clone_samples = defaultdict(list)

# Read in spreadsheet and convert to dict
config_f = glob.glob(f"{args.input_dir}/*.xlsx")[0]
config_f_basename = config_f.rstrip('.xlsx')
sheet = pd.read_excel(config_f, index_col=0)
sheet = sheet.fillna('')
d = sheet.to_dict(orient='index')

clone_for_every_sample = True
for sample_idx in d:
    naming = name(d[sample_idx]['fastq (forward)'])
    sample_dir = os.path.join(args.input_dir, naming)
    pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)
    fq_clone = d[sample_idx]['Clone ID (optional)']

    if fq_clone == '':
        clone_for_every_sample = False
    else:
        clone_samples[fq_clone].append(sample_idx)
        fq_new_base = f"{sample_idx}_{fq_clone}"

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
    reference_path = os.path.join(args.input_dir, 'reference', d[sample_idx]['GOI reference file'])
    # Call pipeline.py for this sample
    subprocess.call(shlex.split(f"python3 {os.path.join(current_dir, 'pipeline.py')} -r {reference_path} -f1 {fq_1_new_path} -f2 {fq_2_new_path}"))
    
    snp_vcf_folder = os.path.join(sample_dir, 'varscan-output')
    snp_vcf_path = os.path.join(snp_vcf_folder, f"{naming}.snp.vcf")
    snp_vcf_aa_path = os.path.join(snp_vcf_folder, f"{naming}.snp.amino_annotated.vcf")
    # Run amino acid annotation on varscan snp output
    subprocess.call(shlex.split(f"python3 {os.path.join(current_dir, '..', 'amino-anno', 'from_varscan.py')} {snp_vcf_path} {reference_path} -o {snp_vcf_folder}"))

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