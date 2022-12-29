
# Cell Line Development NGS pipelines

## Repository overview
This repository contains code for two different variant-calling pipelines: amplicon and RNA-seq. The amplicon pipeline identifies variants and indels from PCR amplicon data, while the RNA-seq library identifies variants and indels from RNA-seq data.

Documentation for each pipeline is in the respective directories: `amp-pip/` for amplicon, and `RNAseq-pip/` for RNA-seq.

Further, there is a third directory, `amino-anno/`, which contains the code for amino acid annotation of VarScan snp reports. This code is used in both the amplicon and RNA-seq pipelines.

### Tagged versions
There are tags to keep track of the work of different interns over the years. `v3-zack` is for work concluded in December 2022 by Zack Mudge. `v2-tongyu-shivam` contains older versions by Tongyu Liu and Shivam Sharma.

## Instructions common between the pipelines
### Environment setup
#### Miniconda setup
There is a miniconda installation in `/mnt/efs/home/cui3008/zack_miniconda/miniconda3`, which contains the `variant-calling` and `RNAseq` environments to run these pipelines. These environments already exist on the Slurm HPC.

If this conda installation isn't currently in your PATH:
1. `cd /mnt/efs/home/cui3008/zack_miniconda/miniconda3/bin`
2. `./conda init bash`. This will configure the conda installation in your `~/.bashrc` file. Close and reopen your ssh connection for the changes to take effect.

### Converting SnapGene .dna file to .fasta
This may be necessary if a reference sequence e.g. GOI sequence is provided as a `.dna` file instead of `.fasta`.

1. Install biopython package with `conda install -c conda-forge biopython`
2. Example Python code for converting `ADT3-2_HC.dna` to `ADT3-2_HC.fasta`:
	```
	from Bio import SeqIO
	records = SeqIO.parse('./ADT3-2_HC.dna', 'snapgene')
	SeqIO.write(records, 'ADT3-2_HC.fasta', 'fasta')
	```
3. Open new `.fasta` file in a text editor and rename the header line accordingly.
