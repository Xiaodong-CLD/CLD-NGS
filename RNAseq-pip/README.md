# Welcome to Cell Line Development RNA-seq pipeline

In this documentation, we will go through steps to:

  1. Set up a working environment.
  2. How to prepare data for your pipeline (refer to amplicon-seq).
  3. Run pipeline. 
  4. Example code.
  5. Generate customized reference genome for future usage. 

This instruction omits the steps of accessing HPC cluster described in the main tutorial. 

### 1. Setting up the working environment.

###### 1.1 Installing the Tools


	conda create -n RNAseq
	conda activate RNAseq
	conda install python=3.7.6
	conda install -c r r
	conda install -c r r-ggplot2
	conda install -c bioconda samtools=1.9
	conda install -c bioconda star=2.7.9
	conda install -c bioconda fastp=0.20.1
	conda install -c bioconda varscan=2.4.4
	conda install -c bioconda pindel=0.2.5b9
	conda install -c bioconda htseq=0.11.3
	conda install -c bioconda picard=2.18.29
	
   ###### NOTE - From this point on, whenever you open a new putty window and log into HPC server, please make sure to type in 
Whenever you are running your pipeline, please make sure you see *(RNAseq)* written at the begining of your username screen.

It should look like: ```(RNAseq) [cui3008@insight-hpc-login2 ~]$```. This *(RNAseq)* at the beginning is extremly important.


If all worked fine, your username line would look like:

```(RNAseq) [cui3008@insight-hpc-login2 ~]$```




    
    
### 2. Prepare Data for your Pipeline.
###### 2.1 Locate pipeline
There are two ways to get pipeline: 

1) Directly execute the scripts saved in ```/mnt/efs/home/cui3008/voof```  

2) Download the scripts from Bitbucket as described in 2.1 (omit in this page)

	
For RNA-seq analysis, the executable script is 
    ```/mnt/efs/home/cui3008/voof/RNAseq-pip/RoarChim.sh``` 


###### 2.3 Organize data files

Our goal is to organize two reads of fastq files from each sample to an individual folder, which is named as sample name. Please refer the tutorial in amplicon-seq or do it manually. 

Your final file name should be named as : *cell-line-name.R1.fastq.gz*, *cell-line-name.R2.fastq.gz*

Also, please make sure that the name: *cell-line-name.R1.fastq.gz* matches the directory name in which it exists.


### 3. Running the Variant Calling Pipeline.

###### 3.1 Standard Running & Testing the Pipeline
This pipeline is designed for gene of interest ADAMTS13. If we want to apply it to other gene of interest, please refer to *Session 4. Generate cutomized reference genome for future usage*. 


		/mnt/efs/home/cui3008/voof/RNAseq-pip/RoarChim.sh \
		-d /path/to/sample/fastq/folder/HD-N62-F23 \
		-g /mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/cho_k1_ncbi_star_2.7 \
		-o /path/to/output/RNA_pipeline/output

 -d is the directory of two fastq files for each sample.
-g is the STAR-built reference genome. You don't have to change it for ADAMTS13 runs.
-o is the user-defined output directory. 
 
A *Portable Bash Script (PBS)* is what you need to access Takeda's HPC. A script has already been written for you.

The script can be submitted using the command:
	```qsub pbs-script.pbs```

### 4. Example code to run pipeline
Path of execuable pipeline.py
 ```/mnt/efs/home/cui3008/voof/RNAseq-pip/RoarChim.sh```
 Path of reference fasta
 ```/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/cho_k1_ncbi_star_2.7```
 Path of the folder stores two fastq files for one sample
 ```/mnt/efs/home/oqh2473/pipeline_test/RNA_pipeline/HD-N62-F23```

Then the command to run the code should be : 

	/mnt/efs/home/cui3008/voof/RNAseq-pip/RoarChim.sh \
	-d /mnt/efs/home/oqh2473/pipeline_test/RNA_pipeline/HD-N62-F23 \
	-g /mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/cho_k1_ncbi_star_2.7 \
	-o /mnt/efs/home/oqh2473/pipeline_test/RNA_pipeline/output
You may not be able to run this code since you don't have the permission to write in this folder. You can check the files in this folder to understand how the intermediate files saved. 

### 5. Generate customized reference genome for future usage.
To generate a customized STAR-built reference genome, the following files are needed (I use ADAMTS13 as example): 
1. The sequence of gene of interest fasta file
```/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/ADAMTS13.fa```
2. The annotation of gene of interest (.gtf). 
```/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/ADAMTS13.gtf```
This file need to be manually write. You can refer ADAMTS13 format. 
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
	


	
