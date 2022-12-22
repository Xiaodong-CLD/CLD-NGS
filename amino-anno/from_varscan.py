#!/usr/bin/env python3

import argparse
import os
import shlex
import subprocess

aa_codes = {
    'A' : 'Ala',
    'C' : 'Cys',
    'D' : 'Asp',
    'E' : 'Glu',
    'F' : 'Phe',
    'G' : 'Gly',
    'H' : 'His',
    'I' : 'Ile',
    'K' : 'Lys',
    'L' : 'Leu',
    'M' : 'Met',
    'N' : 'Asn',
    'P' : 'Pro',
    'Q' : 'Gln',
    'R' : 'Arg',
    'S' : 'Ser',
    'T' : 'Thr',
    'V' : 'Val',
    'W' : 'Trp',
    'Y' : 'Tyr',
    'Z' : 'STOP'
}

# Input: path to fasta file
# Returns: dict with name of sequence from header line, and sequence as one continuous string
def read_fasta(filepath):
    with open(filepath, 'r') as f:
        name = f.readline().rstrip().lstrip('>').split(' ')[0]
        d = {
            'name': name,
            'seq': ''.join(line.rstrip().upper() for line in f)
        }
        return d

# Write sequence to fasta file
# Returns: path of written fasta file
def write_fasta(seq, seqname, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    outpath = os.path.join(out_dir, seqname + '.fa')
    with open(outpath, 'w') as f:
        f.write(f">{seqname}\n")
        f.write(seq + '\n')
    return outpath

# Input: nucleotide sequence as one continuous string
# Returns: amino acid sequence as one continous string
def translate(seq, tmp_dir):
    nt_path = write_fasta(seq, 'nt_tmp', tmp_dir)
    aa_path = os.path.join(tmp_dir, 'aa_tmp.fa')
    subprocess.call(shlex.split(f"{os.path.realpath(os.path.dirname(__file__))}/bin/faTrans {nt_path} {aa_path}"))
    aa_seq = read_fasta(aa_path)['seq']
    os.remove(nt_path)
    os.remove(aa_path)
    return aa_seq

parser = argparse.ArgumentParser()
parser.add_argument('snp', help='Path to varscan snp file')
parser.add_argument('goi', help='Path to GOI FASTA') 
parser.add_argument('-o', help='Output dir', default='./output')
args = parser.parse_args()

if not os.path.exists(args.o):
    os.mkdir(args.o)

ref_fasta = read_fasta(args.goi)
ref_seq = ref_fasta['seq']
ref_trans = translate(ref_seq, args.o)

outpath = os.path.join(args.o, args.snp.split('/')[-1].rsplit('.', 1)[0] + '.amino_annotated.vcf')
with open(args.snp, 'r') as snp_f, open(outpath, 'w') as of:
    fl = snp_f.readline().rstrip().split('\t')
    of.write('\t'.join(fl+['AminoChange', 'MutationType', 'RefAmino', 'VarAmino', 'RefCodon', 'VarCodon']) + '\n')
    for line in snp_f:
        line = line.rstrip().split('\t')
        d = dict(zip(fl, line))

        # Only consider variants in GOI
        ## In RNA-seq pipeline, varscan can identify variants in the host genome. For amino acid annotation, we skip those.
        if d['Chrom'] == ref_fasta['name']:
            idx = int(d['Position']) - 1
            orig_nt = ref_seq[idx]
            seq_l = list(ref_seq)
            seq_l[idx] = d['VarAllele']
            var_seq = ''.join(seq_l)
            var_trans = translate(var_seq, args.o)

            # Determine if amino acid change
            amino_idx = idx // 3
            ref_aa = ref_trans[amino_idx]
            var_aa = var_trans[amino_idx]
            codon_start_idx = idx - idx%3
            ref_codon = ref_seq[codon_start_idx : codon_start_idx+3]
            var_codon = var_seq[codon_start_idx : codon_start_idx+3]
            if ref_aa != var_aa:    # If amino change, determine if nonsense or missense mutation
                d['AminoChange'] = 'true'
                if ref_aa != 'Z' and var_aa == 'Z':
                    d['MutationType'] = 'nonsense'
                elif ref_aa =='Z' and var_aa != 'Z':
                    d['MutationType'] = 'nonstop'
                else:
                    d['MutationType'] = 'missense'
            else:
                d['AminoChange'] = 'false'
                d['MutationType'] = 'N/A'
            d['RefAmino'] = aa_codes[ref_aa]
            d['VarAmino'] = aa_codes[var_aa]
            d['RefCodon'] = ref_codon
            d['VarCodon'] = var_codon
            of.write('\t'.join(list(d.values())) + '\n')
