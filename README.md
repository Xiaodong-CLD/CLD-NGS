
# Welcome to Cell Line Development NGS pipeline

Our current pipeline includes **amplicon-seq pipeline** and **RNA-seq pipeline**.

In this documentation, we will go through steps to:
	1. Set up a working environment.
	2. How to prepare data for your pipeline.
	3. How to run the pipeline.
	4. Possible modifications for the pipeline.


### 1. Setting up the working environment.

Bioinformatics pipelines and applications are best run on Linux environment, manageble on MacOs systems, and not recommended for Windows System. 
Takeda provides access to a High Performance Computing (HPC) cluster to its employees, and an access to such a system is possible through your very own laptop/desktop (which most probably runs on a Windows operating system). 
More information on Takeda's HPC can be found on: [Takeda Insight](https://conf-rds.takeda.com/display/DEET/Insight)

###### 1.1 Accessing Takeda's HPC Cluster/Remote Systems.
In order to access the HPC, you will first need to install the software: *Putty* on your system.
This software lets you connect to other remote (and sometimes more powerful) computers through your laptop.

You can install *Putty* on your system from [Takeda SharePoint](https://mytakeda.sharepoint.com/sites/mytakeda/sitepages/en-us/home.aspx) -> AskIT ticket. Select *Putty 0.70*.

###### 1.2 Getting the right working Environment.
Once you have *Putty* on your system, you need to log-in to Takeda's HPC cluster. At the time of writing this document, the Node Address is: *US HPC LOGIN NODE: 10.245.45.40* [Takeda HPC](https://conf-rds.takeda.com/display/DEET/HPC
You enter ```CUI3008@hpclogin1.rda.onetakeda.com``` in *Host Name* form field and click on *Open*. Enter you password and you should see something like:
		```cui3008@insight-hpc-login2  $```
        
Although, this pipeline can run on almost all the Unix/Linux environment, it is highly recommened to set up a *conda* environment on your Linux system.
In order to get *conda* on your system, go to [Miniconda's website](https://docs.conda.io/en/latest/miniconda.html). Right click on *Miniconda3 Linux 64-bit* for *Python 3* in the *Linux Installers* section; and select *Copy Link Address*.

Go back to you Putty window and write:
	```wget [CTRL+V  or  Right Click to Paste the Link]```
    
You should have successfully downloaded *conda* on your system. Now we will try to install *conda* on you HPC cluster account space. Please write:
	```bash Miniconda-Linux-64...```
  Please follow the instructions of the script and you will have to put in *Yes* in a couple of places.
  
Once you have completed this step, please close the Putty window and open Putty again and log into HPC again.

###### 1.3 Installing the Tools

Once you log back in, you should see something like:

```(base) [cui3008@insight-hpc-login2 ~]$```

Obtain the package file from Lijuan, and run

	conda create -n variant-calling --file package-list.txt

###### NOTE - From this point on, whenever you open a new putty window and log into HPC server, please make sure to type in 

	conda activate variant-calling
    
Whenever you are running your pipeline, please make sure you see *(variant-calling)* written at the begining of your username screen.

It should look like: ```(variant-calling) [cui3008@insight-hpc-login2 ~]$```. This *(variant-calling)* at the beginning is extremly important.


If all worked fine, your username line would look like:

```(variant-calling) [cui3008@insight-hpc-login2 ~]$```


and so on...

After this: Please also run the following command to save coverage JPEG plots.
	
    pip install -U kaleido
    
    
### 2. Prepare Data for your Pipeline.
There are two ways to get pipeline: 
1) Directly execute the scripts saved in ```/mnt/efs/home/cui3008/voof```   
2) Download the scripts from Bitbucket as described in 2.1

###### 2.1 Getting GIT on your system

GIT is a version control system used by programmers all across the world. First make sure that you have GIT on your system by typing: 
	```git```
    
You can go to the following link to visit Takeda's own system to keep track of code repositories.
[Takeda's Bitbucket](https://bb-rds.takeda.com/dashboard)

Our code for pipeline is stored in here, and in order to access this, you will need to set up an SSH key on your system (SSH key is like a password that lets Takeda's Bitbucket to recognize you). 
We'll follow this [Tutorial](https://docs.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent) little tutorial by GitHub to get this SSH key.

Please type:

	ssh-keygen -t rsa -b 4096 -C "your_email@example.com"

on your system and don't forget to put your Takeda email in the last part of the above command. Type in "Yes" wherever necessary.

If asked for *Passphrase*, you can enter a simple password (don't forget it), or even leave it blank. It is just an additional security step.

Once this is done, type in:

	cat ~/.ssh/id_rsa.pub
    
You'll see around 3-5 lines of random characters (text). Please select this whole text and copy it using ctrl+c. The copied text should start from ssh-rsa and end at your.email@takeda.com


Now you need to go to [Takeda's Bitbucket](https://bb-rds.takeda.com/dashboard), and click on top right user profile icon to select *Manage Account*.

Please select *SSH Keys* from the left column and click on *Add Key* and paste your copied key text there.
    

###### 2.2 Installing the Tools
The second step is to get the code base for you pipeline. We will use *GIT* software to manage everything on the system. I would expect *GIT* to be already there on your system, but if it isn't there, then you can always install it using:

	conda install -c anaconda git
    
Then, to get the code repository:

	git clone ssh://git@bb-rds.takeda.com:8080/cldngs/voof.git
    
If you get a permission error, then please ask Dr. Shen to give you a permission to get the repo. But if you are able to access this document, you most probably already have the permission to run the command.

You can check everything was alright using the command:
	```ls```

Use the following command to enter the code repository:
	```cd voof```
    
The key directory of interest is *pipeline*, so go to pipeline directory using:

For amplicon-seq
	```cd amp-pip/``` 
	
For RNA-seq analysis
    ```cd RNAseq-pip/``` 
###### 2.3 Run the tool requirement check file here.

In order to run the software requirement checker, type in the following command:
 	```python3 check_pipeline_tools_requirement.py```

If the last output line is: *[Takeda Pipeline Message] Go AHEAD! All the required tools are on the system.*,
you're all set! If you see a tool missing, please install it using conda. A quick google search will give you the exact command for installation of any bioinformatics tool.


###### 2.4 Create the Data Repository. 

Make directory for your reference files and fastq files by ```mkdir /path/to/fastq```. You can use ```pwd``` to check the full path of your working directory. 

Cyberduck is a good file transfer tool to copy and paste data from your local disk drive to HPC. You can submit a ticket to AskIT to install it. 

*Note:* We usually do not store data on Bitbucket/GIT, and so you'll have to transfer the data using other means. We're still working on managing optimal sharing of data.


###### 2.5 Organize data files

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
	cp $c* /path/to/Fastq/$c
	mv *R1*gz $c.R1.fastq.gz
	mv *R2*gz $c.R2.fastq.gz
	done

If you are not comfortable with this method, you can always do it manually (please refer the exact steps in previous pipeline tutorial). 

Your final file name should be named as : *cell-line-name.R1.fastq.gz*, *cell-line-name.R2.fastq.gz*
Also, please make sure that the name: *cell-line-name.R1.fastq.gz* matches the directory name in which it exists.


### 3. Running the Variant Calling Pipeline.

###### 3.1 Standard Running & Testing the Pipeline

There are several aspects and ways in which you can run the pipeline. But the easiest way is using the following command:

```python3 /mnt/efs/home/cui3008/voof/amp-pip/pipeline.py -r /path/to/reference.fasta -f1 /path/to/sample-name.R1.fastq.gz -f2 /path/to/sample-name.R1.fastq.gz```

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
