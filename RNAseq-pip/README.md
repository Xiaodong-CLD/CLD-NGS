# Welcome to Cell Line Development RNA-seq pipeline

In this documentation, we will go through steps to:

  1. Set up a working environment.
  2. How to prepare data for your pipeline.
  3. Run pipeline. 
  4. Generate customized reference genome for future usage. 


### 1. Setting up the working environment.
1. Setup conda installation following the instructions in the main repository README.
2. `conda activate RNAseq`. The `RNAseq` environment already exists on the Slurm HPC, but `RNAseq.yml` contains an export of the packages in case conda needs to be reinstalled.

### 2. Prepare Data for your Pipeline.
###### 2.1 Locate pipeline
For RNA-seq analysis, the executable script is `RoarChim.sh`

###### 2.2 Organize data files

Our goal is to organize two reads of fastq files from each sample to an individual folder, which is named as sample name.

Your final file name should be named as : *cell-line-name.R1.fastq.gz*, *cell-line-name.R2.fastq.gz*

Also, please make sure that the name: *cell-line-name.R1.fastq.gz* matches the directory name in which it exists.


### 3. Running the Variant Calling Pipeline.

###### 3.1 Standard Running & Testing the Pipeline
This pipeline is designed for gene of interest ADAMTS13. If we want to apply it to other gene of interest, please refer to *Section 4. Generate cutomized reference genome for future usage*. 

See `job.script` for an example Slurm script containing the call to `RoarChim.sh`. Submit a job to the queue with `sbatch job.script`.

RoarChim.sh takes three arguments:

	-d is the directory of two fastq files for each sample.
	-g is the STAR-built reference genome. You don't have to change it for ADAMTS13 runs.
	-o is the user-defined directory path to save final outputs. Create this empty directory before running the pipeline.
	
### 4. Generate customized reference genome for future usage.
To generate a customized STAR-built reference genome, the following files are needed (I use ADAMTS13 as example): 

1. The sequence of gene of interest fasta file

```/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/ADAMTS13.fa```

2. The annotation of gene of interest (.gtf). 

```/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/ADAMTS13.gtf```

This file needs to be manually written. You can use ADAMTS13.gtf as template. 

3. The sequence of CHO genome

```/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna``` 

4. The annotation of CHO genome

 ```/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf```
 
CHO genome sequence and annotation are downloaded from 
 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045/GCF_003668045.3_CriGri-PICRH-1.0/

Then we will make customized CHO+goi genome: 

	cat ADAMTS13.fa GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna > adam_cho.fa
	cat ADAMTS13.gtf GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf > adam_cho.gtf
    
Use STAR to build reference genome for alignment:

	STAR --runMode genomeGenerate --runThreadN 8 \
	--genomeDir /path/to/reference/output/cho_k1_ncbi_star \
	--genomeFastaFiles /path/to/adam_cho.fa \
	--sjdbGTFfile /path/to/adam_cho.gtf \
	--sjdbOverhang 149 
--sjdbOverhang was set as readlength -1 

Next, we need to edit RoarChim.sh to change some variables. Please save RoarChim.sh a copy in case we mess up during modification. A bed file will be needed as addition. Please refer the ADAMTS13.bed to make the new one for gene of interest. 

	gtf=/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/ADAMTS13.gtf
	ref_seq=/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/adam_cho.fa
	goi_seq=/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/ADAMTS13.fa
	bed=/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/ADAMTS13.bed

Change the above variables as needed.

Then scan the script to change anywhere of 'ADAMTS13' to the name of gene of interest. 

At last, run the pipeline with the new RoarChim.sh with your fastq directory and new reference genome path.  

For any questions, feel free to email Ltongyu@umich.edu
	


	
