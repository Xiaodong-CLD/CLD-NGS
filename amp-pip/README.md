
# Welcome to Cell Line Development amplicon-seq pipeline

In this documentation, we will go through steps to:
	1. Set up a working environment.
	2. How to prepare data for your pipeline.
	3. How to run the pipeline.
	4. Example code
This instruction omits the steps of accessing HPC cluster described in the main tutorial. 

### 1. Setting up the working environment.

###### 1.1 Installing the Tools


	conda create -n variant-calling
	conda activate variant-calling
	conda install python=3.7.6
	conda install -c bioconda samtools=1.9
	conda install -c bioconda bwa=0.7.17
	conda install -c bioconda fastp=0.20.1
	conda install -c bioconda varscan=2.4.4
	conda install -c bioconda pindel=0.2.5b9
	conda install -c bioconda bamtools=2.5.1
	
   ###### NOTE - From this point on, whenever you open a new putty window and log into HPC server, please make sure to type in 
Whenever you are running your pipeline, please make sure you see *(variant-calling)* written at the begining of your username screen.

It should look like: ```(variant-calling) [cui3008@insight-hpc-login2 ~]$```. This *(variant-calling)* at the beginning is extremly important.


If all worked fine, your username line would look like:

```(variant-calling) [cui3008@insight-hpc-login2 ~]$```


and so on...

After this: Please also run the following command to save coverage JPEG plots.
	
	pip install plotly
	pip install -U kaleido

    
    
### 2. Prepare Data for your Pipeline.
###### 2.1 Locate pipeline
There are two ways to get pipeline: 
1) Directly execute the scripts saved in ```/mnt/efs/home/cui3008/voof```   
2) Download the scripts from Bitbucket as described in 2.1 (omit in this page)

Use the following command to enter the code repository:
	```cd /mnt/efs/home/cui3008/voof```
    
The key directory of interest is *pipeline*, so go to pipeline directory using:

For amplicon-seq
	```cd ./amp-pip/``` 
	
For RNA-seq analysis
    ```cd ./RNAseq-pip/``` 
###### 2.2 Run the tool requirement check file here.

In order to run the software requirement checker, type in the following command:
 	```python3 check_pipeline_tools_requirement.py```

If the last output line is: *[Takeda Pipeline Message] Go AHEAD! All the required tools are on the system.*,
you're all set! If you see a tool missing, please install it using conda. A quick google search will give you the exact command for installation of any bioinformatics tool.


###### 2.3 Organize data files

Our goal is to organize two reads of fastq files from each sample to an individual folder, which is named as sample name. 

A little magic you can use here: 

Step 1: Transfer the Fastq folder from scientific drive to HPC. 

```cd /path/to/fastq```
By ```ls``` you will see a list of fastq files in the folder

Step 2: Use the code below to move the two fastq files. You need to change ```/path/to/Fastq/``` to your Fastq folder path at first. 

	for f in /path/to/Fastq/*; do
    b=`basename $f`
    echo ${b%%\_*}
    c=${b%%\_*}
    mkdir /path/to/Fastq/$c
    mv $c* /path/to/Fastq/$c
    cd /path/to/Fastq/$c
    mv *R1*gz $c.R1.fastq.gz
    mv *R2*gz $c.R2.fastq.gz
    done

If you are not comfortable with this method, you can always do it manually (please refer the exact steps in previous pipeline tutorial). 

Your final file name should be named as : *cell-line-name.R1.fastq.gz*, *cell-line-name.R2.fastq.gz*
Also, please make sure that the name: *cell-line-name.R1.fastq.gz* matches the directory name in which it exists.


### 3. Running the Variant Calling Pipeline.

###### 3.1 Standard Running & Testing the Pipeline

There are several aspects and ways in which you can run the pipeline. But the easiest way is using the following command:

```python3 /mnt/efs/home/cui3008/voof/amp-pip/pipeline.py -r /path/to/reference.fasta -f1 /path/to/sample-name.R1.fastq.gz -f2 /path/to/sample-name.R2.fastq.gz```

This will show a lot of text/messages on your screen and can take a lot of time to run depending upon the system you have.

###### 3.2 Using the pipeline_runner.py Script

###### 3.3 Running the application on the HPC cluster using PBS script.

A *Portable Bash Script (PBS)* is what you need to access Takeda's HPC. A script has already been written for you.
The script can be submitted using the command:
	```qsub pbs-script.pbs```

### 4. Example code to run pipeline
Path of execuable pipeline.py
 ```/mnt/efs/home/cui3008/voof/amp-pip/pipeline.py```
 Path of reference fasta
 ```/mnt/efs/home/oqh2473/pipeline_test/VWFcp_ADAMTS13.fasta```
 Path of two fastq files
 ```/mnt/efs/home/cui3008/amp_202108/analysis/F6-1/F6-1.R1.fastq.gz```
 ```/mnt/efs/home/cui3008/amp_202108/analysis/F6-1/F6-1.R2.fastq.gz```
 
Then the command to run the code should be : 

	python3 /mnt/efs/home/cui3008/voof/amp-pip/pipeline.py \
	-r  /mnt/efs/home/oqh2473/pipeline_test/VWFcp_ADAMTS13.fasta \
	-f1 /mnt/efs/home/cui3008/amp_202108/analysis/F6-1/F6-1.R1.fastq.gz \
	-f2 /mnt/efs/home/cui3008/amp_202108/analysis/F6-1/F6-1.R2.fastq.gz
You may not be able to run this code since you dont have the permission to write in this folder. You can check the files in this folder to understand how the intermediate files saved. 



