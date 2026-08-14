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
# Also writes a run-level multi-caller indel summary (Results/indel.txt)
# combining Pindel, VarScan, and LoFreq calls; only positions supported by
# ≥2 callers with ≥1% read support are reported.

def name(string):
    prefix = string.split("_L")[0]
    return prefix


# Minimum AD fractional support for Pindel calls: alt / (ref + alt) >= 1%.
# Pindel AD overstates true VAF ~4×, so this is not 1% true AF. See Exp C.
PINDEL_MIN_AD_FRAC = 0.01


def parse_pindel_summary(summary_path, sample_id, clone_id):
    # Parse a per-sample Pindel Summary.vcf into a tidy DataFrame so the
    # structural-variant (indel) calls can be merged into the run-level
    # Results/indel.txt report, laid out like snp.txt.
    # Only rows with AD alt / (ref + alt) >= PINDEL_MIN_AD_FRAC are kept.
    rows = []
    if not os.path.exists(summary_path):
        return pd.DataFrame(rows)
    with open(summary_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            f = line.split('\t')
            if len(f) < 8:
                continue
            info = dict(
                kv.split('=', 1) if '=' in kv else (kv, '')
                for kv in f[7].split(';')
            )
            gt = ref_reads = alt_reads = var_freq = ''
            ad_frac = 0.0
            if len(f) >= 10 and ':' in f[8]:
                fd = dict(zip(f[8].split(':'), f[9].split(':')))
                gt = fd.get('GT', '')
                if ',' in fd.get('AD', ''):
                    ref_reads, alt_reads = fd['AD'].split(',')[:2]
                    try:
                        r, a = int(ref_reads), int(alt_reads)
                        if r + a > 0:
                            ad_frac = a / (r + a)
                            var_freq = f"{round(100 * ad_frac, 2)}%"
                    except ValueError:
                        pass
            if ad_frac < PINDEL_MIN_AD_FRAC:
                continue
            rows.append({
                'SampleID': sample_id,
                'CloneID': clone_id,
                'Chrom': f[0],
                'Position': f[1],
                'Ref': f[3],
                'Alt': f[4],
                'SVTYPE': info.get('SVTYPE', ''),
                'SVLEN': info.get('SVLEN', ''),
                'HOMLEN': info.get('HOMLEN', ''),
                'GT': gt,
                'RefReads': ref_reads,
                'AltReads': alt_reads,
                'VarFreq_pindel': var_freq,
            })
    return pd.DataFrame(rows)


def parse_varscan_indels(indel_vcf_path, sample_id, clone_id):
    # Parse VarScan pileup2indel TSV output into a tidy DataFrame.
    # VarScan indel TSV columns (0-indexed): Chrom, Position, Ref, Cons,
    # Reads1, Reads2, VarFreq, ... (same layout as pileup2snp, 15 cols).
    # Only rows with VarFreq >= 1% are kept.
    rows = []
    if not os.path.exists(indel_vcf_path):
        return pd.DataFrame(rows)
    with open(indel_vcf_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            f = line.split('\t')
            if len(f) < 7:
                continue
            # Skip the header row (Chrom column header is literally "Chrom").
            if f[0] == 'Chrom':
                continue
            var_freq_str = f[6]  # e.g. "12.34%"
            try:
                var_freq = float(var_freq_str.rstrip('%'))
            except ValueError:
                continue
            if var_freq < 1.0:
                continue
            rows.append({
                'SampleID': sample_id,
                'CloneID': clone_id,
                'Chrom': f[0],
                'Position': f[1],
                'Ref': f[2],
                'Alt': f[3],  # Cons, e.g. "+1T" or "-3"
                'VarFreq_varscan': f"{round(var_freq, 2)}%",
            })
    return pd.DataFrame(rows)


def parse_lofreq_indels(lofreq_vcf_path, sample_id, clone_id):
    # Parse LoFreq VCF indel output into a tidy DataFrame.
    # LoFreq INFO field contains AF=<float> (allele frequency, 0-1).
    # Only true indel rows (len(REF) != len(ALT)) with AF >= 1% are kept.
    rows = []
    if not os.path.exists(lofreq_vcf_path):
        return pd.DataFrame(rows)
    with open(lofreq_vcf_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            f = line.split('\t')
            if len(f) < 8:
                continue
            chrom, pos, ref, alt = f[0], f[1], f[3], f[4]
            # Confirm this is an indel (lengths differ).
            if len(ref) == len(alt):
                continue
            info = dict(
                kv.split('=', 1) if '=' in kv else (kv, '')
                for kv in f[7].split(';')
            )
            try:
                var_freq = round(float(info.get('AF', 0)) * 100, 2)
            except ValueError:
                continue
            if var_freq < 1.0:
                continue
            rows.append({
                'SampleID': sample_id,
                'CloneID': clone_id,
                'Chrom': chrom,
                'Position': pos,
                'Ref': ref,
                'Alt': alt,
                'VarFreq_lofreq': f"{var_freq}%",
            })
    return pd.DataFrame(rows)


def merge_indel_callers(pindel_df, varscan_df, lofreq_df):
    # Merge per-sample indel calls from Pindel, VarScan, and LoFreq.
    # Merge key: (SampleID, Chrom, Position).
    # Returns all rows detected by any caller at ≥1% AF, with a 'source'
    # column listing which caller(s) found it (e.g. "pindel", "pindel,lofreq").

    # Build a combined record dict: key -> {caller freq fields + Ref/Alt/CloneID}.
    combined = {}
    for caller, df, freq_col in [
        ('pindel', pindel_df, 'VarFreq_pindel'),
        ('varscan', varscan_df, 'VarFreq_varscan'),
        ('lofreq', lofreq_df, 'VarFreq_lofreq'),
    ]:
        if df.empty:
            continue
        for _, row in df.iterrows():
            key = (row['SampleID'], row['Chrom'], row['Position'])
            if key not in combined:
                combined[key] = {
                    'SampleID': row['SampleID'],
                    'CloneID': row.get('CloneID', ''),
                    'Chrom': row['Chrom'],
                    'Position': row['Position'],
                    'Ref': row.get('Ref', ''),
                    'Alt': row.get('Alt', ''),
                }
            combined[key][freq_col] = row[freq_col]

    rows = []
    for data in combined.values():
        callers = [c for c in ['pindel', 'varscan', 'lofreq']
                   if f'VarFreq_{c}' in data]
        data['source'] = ','.join(callers)
        rows.append(data)

    return pd.DataFrame(rows)


current_dir = pathlib.Path(__file__).parent.resolve()

parser = argparse.ArgumentParser()
parser.add_argument('input_dir', help='Path to input directory')
args = parser.parse_args()

snp_dfs = []
pindel_dfs = []
varscan_indel_dfs = []
lofreq_indel_dfs = []
clone_samples = defaultdict(list)

# Read in spreadsheet and convert to dict
config_f = glob.glob(f"{args.input_dir}/sample_xlsx/*.xlsx")[0]
config_f_basename = config_f.rstrip('.xlsx')
sheet = pd.read_excel(config_f, index_col=0)
sheet = sheet.fillna('')
d = sheet.to_dict(orient='index')

folders_list = list()
clone_for_every_sample = True
for sample_idx in d:
    naming = name(d[sample_idx]['fastq (forward)'])
    # Results/<sample>/ — user-facing outputs only (VCFs, plots, QC reports).
    sample_dir = os.path.join(args.input_dir, "Results", naming)
    # tmp/<sample>/ — all pipeline intermediates (BAMs, mpileup, raw VCFs, pindel intermediates).
    tmp_dir = os.path.join(args.input_dir, "tmp", naming)
    pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)
    pathlib.Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    folders_list.append(sample_dir)

    fq_clone = d[sample_idx].get('Clone ID (optional)', '')

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
    # Symlinks go to tmp_dir so all derived files (trimmed FASTQs, BAMs, etc.) land there.
    fq_1_new_path = os.path.join(tmp_dir, fq_1_new_name)
    fq_2_new_path = os.path.join(tmp_dir, fq_2_new_name)
    # Symlink the renamed FASTQs instead of copying (avoids a per-run full-file copy).
    for src, dst in ((fq_1_path, fq_1_new_path), (fq_2_path, fq_2_new_path)):
        if os.path.islink(dst) or os.path.exists(dst):
            os.remove(dst)
        os.symlink(os.path.abspath(src), dst)
    reference_path = os.path.join(args.input_dir, 'reference', d[sample_idx]['GOI reference file'])
    # Call pipeline.py for this sample
    subprocess.call(shlex.split(f"python3 {os.path.join(current_dir, 'pipeline.py')} -r {reference_path} -f1 {fq_1_new_path} -f2 {fq_2_new_path}"))
    
    # Amino annotation reads from tmp_dir where pipeline.py wrote the raw SNP VCF.
    snp_vcf_folder = os.path.join(tmp_dir, 'varscan-output')
    snp_vcf_path = os.path.join(snp_vcf_folder, f"{naming}.snp.vcf")
    snp_vcf_aa_path = os.path.join(snp_vcf_folder, f"{naming}.snp.amino_annotated.vcf")
    # Run amino acid annotation on varscan snp output
    subprocess.call(shlex.split(f"python3 {os.path.join(current_dir, '..', 'amino-anno', 'from_varscan.py')} {snp_vcf_path} {reference_path} -o {snp_vcf_folder}"))

    # ── Copy user-facing outputs from tmp_dir to Results/<sample>/ ──────────
    # fastp QC reports
    for p in (glob.glob(os.path.join(tmp_dir, "*.fastp.html"))
              + glob.glob(os.path.join(tmp_dir, "*.fastp.json"))):
        shutil.copy2(p, sample_dir)
    # Coverage plots
    cov_src = os.path.join(tmp_dir, "coverage-plot-output")
    if os.path.isdir(cov_src):
        shutil.copytree(cov_src, os.path.join(sample_dir, "coverage-plot-output"),
                        dirs_exist_ok=True)
        # Also copy into run-level aggregated folder (Results/coverage-plot-output/).
        agg_cov_dir = os.path.join(args.input_dir, "Results", "coverage-plot-output")
        pathlib.Path(agg_cov_dir).mkdir(exist_ok=True)
        shutil.copytree(cov_src, agg_cov_dir, dirs_exist_ok=True)
    # VarScan SNP (amino-annotated) and indel VCFs
    results_varscan_dir = os.path.join(sample_dir, "varscan-output")
    pathlib.Path(results_varscan_dir).mkdir(exist_ok=True)
    for fname in (f"{naming}.snp.vcf", f"{naming}.indel.vcf",
                  f"{naming}.snp.amino_annotated.vcf"):
        src = os.path.join(snp_vcf_folder, fname)
        if os.path.exists(src):
            shutil.copy2(src, results_varscan_dir)
    # Pindel final Summary VCF
    results_pindel_dir = os.path.join(sample_dir, "pindel-output")
    pathlib.Path(results_pindel_dir).mkdir(exist_ok=True)
    pindel_summary_src = os.path.join(tmp_dir, "pindel-output", "pindel-output.Summary.vcf")
    if os.path.exists(pindel_summary_src):
        shutil.copy2(pindel_summary_src, results_pindel_dir)
    # LoFreq indel VCF
    results_lofreq_dir = os.path.join(sample_dir, "lofreq-output")
    pathlib.Path(results_lofreq_dir).mkdir(exist_ok=True)
    lofreq_vcf_src = os.path.join(tmp_dir, "lofreq-output", f"{naming}.lofreq.indel.vcf")
    if os.path.exists(lofreq_vcf_src):
        shutil.copy2(lofreq_vcf_src, results_lofreq_dir)
    # ────────────────────────────────────────────────────────────────────────

    snp_df = pd.read_csv(snp_vcf_aa_path, sep='\t')
    snp_df.insert(0, 'SampleID', sample_idx)
    snp_df['CloneID'] = fq_clone
    # Keeping snp output dict for each sample
    snp_dfs.append(snp_df)

    # Collect per-sample indel calls from all three callers.
    pindel_summary_path = os.path.join(tmp_dir, 'pindel-output', 'pindel-output.Summary.vcf')
    pindel_dfs.append(parse_pindel_summary(pindel_summary_path, sample_idx, fq_clone))
    varscan_indel_path = os.path.join(tmp_dir, 'varscan-output', f'{naming}.indel.vcf')
    varscan_indel_dfs.append(parse_varscan_indels(varscan_indel_path, sample_idx, fq_clone))
    lofreq_indel_path = os.path.join(tmp_dir, 'lofreq-output', f'{naming}.lofreq.indel.vcf')
    lofreq_indel_dfs.append(parse_lofreq_indels(lofreq_indel_path, sample_idx, fq_clone))

# Putting all the individual snp dicts together
snp_combined_orig = pd.concat(snp_dfs)
# Dropping unnecessary columns
snp_combined_orig = snp_combined_orig.drop(['Strands1', 'Strands2', 'Qual1', 'Qual2', 'MapQual1', 'MapQual2'], axis=1)
# Converting string percentages to float
snp_combined_orig['VarFreq'] = snp_combined_orig['VarFreq'].apply(lambda x: float(x[:-1]))
# Dropping all entries with VarFreq < 1%
snp_combined_orig = snp_combined_orig[snp_combined_orig['VarFreq'] >= 1]
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

    with open(f"{args.input_dir}/Results/snp.txt", 'w') as out_f:
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
                out_f.write('No positions with variants above 1% frequency in all replicates\n\n')
            else:
                out_f.write(f"{clone_output_df.to_string(index=False)}\n\n")
else:
    with open(f"{args.input_dir}/Results/snp.txt", 'w') as out_f:
        for sample_id in d:
            out_f.write(f"SampleID: {sample_id}\n")
            sample_df = snp_combined_orig.loc[snp_combined_orig['SampleID'] == sample_id] 
            if sample_df.empty:
                out_f.write('No positions with variants above 1% frequency\n\n')
            else:
                out_f.write(f"{sample_df.to_string(index=False)}\n\n")

# Merge per-sample indel calls from Pindel, VarScan, and LoFreq.
# Results/indel.txt reports all positions detected by any caller at ≥1% support.
pindel_all = pd.concat(pindel_dfs, ignore_index=True) if pindel_dfs else pd.DataFrame()
varscan_all = pd.concat(varscan_indel_dfs, ignore_index=True) if varscan_indel_dfs else pd.DataFrame()
lofreq_all = pd.concat(lofreq_indel_dfs, ignore_index=True) if lofreq_indel_dfs else pd.DataFrame()
indel_combined = merge_indel_callers(pindel_all, varscan_all, lofreq_all)
indel_display_cols = ['Chrom', 'Position', 'Ref', 'Alt',
                      'VarFreq_pindel', 'VarFreq_varscan', 'VarFreq_lofreq', 'source']
# Replace NaN caller-freq cells with '-' for readability.
if not indel_combined.empty:
    indel_combined = indel_combined.fillna({c: '-' for c in
                                            ['VarFreq_pindel', 'VarFreq_varscan', 'VarFreq_lofreq']})
with open(f"{args.input_dir}/Results/indel.txt", 'w') as out_f:
    if clone_for_every_sample:
        for clone_id, sample_ids in clone_samples.items():
            out_f.write(f"CloneID: {clone_id}\n")
            out_f.write(f"Replicate Sample IDs: {', '.join(str(x) for x in sample_ids)}\n")
            if indel_combined.empty:
                clone_df = indel_combined
            else:
                clone_df = indel_combined[indel_combined['SampleID'].isin(sample_ids)]
            if clone_df.empty:
                out_f.write('No indels detected at ≥1% frequency\n\n')
            else:
                out_f.write(f"{clone_df[['SampleID'] + indel_display_cols].to_string(index=False)}\n\n")
    else:
        for sample_id in d:
            out_f.write(f"SampleID: {sample_id}\n")
            if indel_combined.empty:
                sample_df = indel_combined
            else:
                sample_df = indel_combined[indel_combined['SampleID'] == sample_id]
            if sample_df.empty:
                out_f.write('No indels detected at ≥1% frequency\n\n')
            else:
                out_f.write(f"{sample_df[indel_display_cols].to_string(index=False)}\n\n")

print(f"[Takeda Pipeline Message] Pipeline complete. Results written to {os.path.join(args.input_dir, 'Results')} "
      f"(per-sample folders, snp.txt, indel.txt).")