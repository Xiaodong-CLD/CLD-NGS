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
    sample_dir = os.path.join(args.input_dir, str(sample_idx))
    pathlib.Path.mkdir(sample_dir)
    
    fq_clone = d[sample_idx]['Clone ID (optional)']
    fq_new_base = str(sample_idx)
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
    fq_1_new_name = f"{fq_new_base}.R1{ext_1}{fq_1_split[1]}"

    fq_2_name = d[sample_idx]['fastq (reverse)']
    if '.fastq' in fq_2_name:
        ext_2 = '.fastq'
    else:
        ext_2 = '.fq'
    fq_2_split = fq_2_name.split(ext_1)
    fq_2_new_name = f"{fq_new_base}.R2{ext_2}{fq_2_split[1]}"

    fq_1_path = os.path.join(args.input_dir, 'fastq', fq_1_name)
    fq_2_path = os.path.join(args.input_dir, 'fastq', fq_2_name)
    fq_1_new_path = os.path.join(sample_dir, fq_1_new_name)
    fq_2_new_path = os.path.join(sample_dir, fq_2_new_name)
    shutil.copyfile(fq_1_path, fq_1_new_path)
    shutil.copyfile(fq_2_path, fq_2_new_path)
    reference_path = os.path.join(args.input_dir, 'reference', d[sample_idx]['GOI reference file'])
    # Call pipeline.py for this sample
    subprocess.call(shlex.split(f"python3 {os.path.join(current_dir, 'pipeline.py')} -r {reference_path} -f1 {fq_1_new_path} -f2 {fq_2_new_path}"))

    snp_vcf_folder = os.path.join(sample_dir, 'varscan-output')
    snp_vcf_path = os.path.join(snp_vcf_folder, f"{fq_new_base}.snp.vcf")
    snp_vcf_aa_path = os.path.join(snp_vcf_folder, f"{fq_new_base}.snp.amino_annotated.vcf")
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
snp_combined = snp_combined_orig.set_index(['SampleID', 'Position'])
snp_combined_orig = snp_combined_orig.drop('CloneID', axis=1)

if clone_for_every_sample:
    # If clone information is provided for every sample,
    # write out snp report considering the replicates
    by_sample_pos = snp_combined.to_dict('index')
    by_clone_by_pos = defaultdict(lambda: defaultdict(list))
    for k,v in by_sample_pos.items():
        by_clone_by_pos[v['CloneID']][k[1]].append(k[0])

    with open(f"{config_f_basename}.snp.txt", 'w') as out_f:
        for clone_id,sample_ids in clone_samples.items():
            clone_output_df = pd.DataFrame()
            sample_l = clone_samples[clone_id]
            pos_samples = by_clone_by_pos[clone_id] 
            for pos in pos_samples:
                if pos_samples[pos] == sample_l:
                    for s in sample_l:
                        row = snp_combined_orig.query("SampleID == @s and Position == @pos")
                        clone_output_df = pd.concat([clone_output_df, row]) 
            out_f.write(f"CloneID: {clone_id}\n")
            out_f.write(f"Replicate Sample IDs: {', '.join(str(x) for x in sample_l)}\n")
            if clone_output_df.empty:
                out_f.write('No positions with variants above 2% frequency in all replicates\n\n')
            else:
                out_f.write(f"{clone_output_df.to_string(index=False)}\n\n")
else:
    # Since clone information wasn't provided for every sample,
    # we write out the report by individual samples
    with open(f"{config_f_basename}.snp.txt", 'w') as out_f:
        for sample_id in d:
            out_f.write(f"SampleID: {sample_id}\n")
            sample_df = snp_combined_orig.loc[snp_combined_orig['SampleID'] == sample_id] 
            if sample_df.empty:
                out_f.write('No positions with variants above 2% frequency\n\n')
            else:
                out_f.write(f"{sample_df.to_string(index=False)}\n\n")
