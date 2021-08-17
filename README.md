Hello, Welcome to Cell Line Development's Variant Calling Pipeline.
In this documentation, we will go through steps to: 1. Set up a working environment. 2. How to prepare data for your pipeline. 3. How to run the pipeline. 4. Possible modifications for the pipeline.

1. Setting up the working environment.
Bioinformatics pipelines and applications are best run on Linux environment, manageble on MacOs systems, and not recommended for Windows System. Takeda provides access to a High Performance Computing (HPC) cluster to its employees, and an access to such a system is possible through your very own laptop/desktop (which most probably runs on a Windows operating system). More information on Takeda's HPC can be found on: Takeda Insight

1.1 Accessing Takeda's HPC Cluster/Remote Systems.
In order to access the HPC, you will first need to install the software: Putty on your system. This software lets you connect to other remote (and sometimes more powerful) computers through your laptop.

You can install Putty on your system from Takeda SharePoint -> AskIT ticket. Select Putty 0.70.

1.2 Getting the right working Environment.
Once you have Putty on your system, you need to log-in to Takeda's HPC cluster. At the time of writing this document, the Node Address is: US HPC LOGIN NODE: 10.245.45.40 [Takeda HPC](https://conf-rds.takeda.com/display/DEET/HPC You enter CUI3008@10.245.45.40 in Host Name form field and click on Open. Enter you password and you should see something like: cui3008@insight-hpc-login2 $

Although, this pipeline can run on almost all the Unix/Linux environment, it is highly recommened to set up a conda environment on your Linux system. In order to get conda on your system, go to Miniconda's website. Right click on Miniconda3 Linux 64-bit for Python 3 in the Linux Installers section; and select Copy Link Address.

Go back to you Putty window and write: wget [CTRL+V or Right Click to Paste the Link]

You should have successfully downloaded conda on your system. Now we will try to install conda on you HPC cluster account space. Please write: bash Miniconda-Linux-64... Please follow the instructions of the script and you will have to put in Yes in a couple of places.

Once you have completed this step, please close the Putty window and open Putty again and log into HPC again.

1.3 Installing the Tools
Once you log back in, you should see something like:

(base) [cui3008@insight-hpc-login2 ~]$

Obtain the package file from Lijuan, and run

conda create -n variant-calling --file package-list.txt
NOTE - From this point on, whenever you open a new putty window and log into HPC server, please make sure to type in
conda activate variant-calling
Whenever you are running your pipeline, please make sure you see (variant-calling) written at the begining of your username screen.

It should look like: (variant-calling) [cui3008@insight-hpc-login2 ~]$. This (variant-calling) at the beginning is extremly important.

If all worked fine, your username line would look like:

(variant-calling) [cui3008@insight-hpc-login2 ~]$

and so on...

After this: Please also run the following command to save coverage JPEG plots.

pip install -U kaleido
2. Prepare Data for your Pipeline.
2.1 Getting GIT on your system
GIT is a version control system used by programmers all across the world. First make sure that you have GIT on your system by typing: git

You can go to the following link to visit Takeda's own system to keep track of code repositories. Takeda's Bitbucket

Our code for pipeline is stored in here, and in order to access this, you will need to set up an SSH key on your system (SSH key is like a password that lets Takeda's Bitbucket to recognize you). We'll follow this Tutorial little tutorial by GitHub to get this SSH key.

Please type:

ssh-keygen -t rsa -b 4096 -C "your_email@example.com"
on your system and don't forget to put your Takeda email in the last part of the above command. Type in "Yes" wherever necessary.

If asked for Passphrase, you can enter a simple password (don't forget it), or even leave it blank. It is just an additional security step.

Once this is done, type in:

cat ~/.ssh/id_rsa.pub
You'll see around 3-5 lines of random characters (text). Please select this whole text and copy it using ctrl+c. The copied text should start from ssh-rsa and end at your.email@takeda.com

Now you need to go to Takeda's Bitbucket, and click on top right user profile icon to select Manage Account.

Please select SSH Keys from the left column and click on Add Key and paste your copied key text there.

2.2 Installing the Tools
The second step is to get the code base for you pipeline. We will use GIT software to manage everything on the system. I would expect GIT to be already there on your system, but if it isn't there, then you can always install it using: conda install -c anaconda git

Then, to get the code repository: git clone ssh://git@bb-rds.takeda.com:8080/var/variant-calling-pipeline.git

If you get a permission error, then please ask Dr. Shen to give you a permission to get the repo. But if you are able to access this document, you most probably already have the permission to run the command.

You can check everything was alright using the command: ls

Use the following command to enter the code repository: cd variant-calling-pipeline

The key directory of interest is pipeline, so go to pipeline directory using: cd pipeline/

2.3 Run the tool requirement check file here.
In order to run the software requirement checker, type in the following command: python3 check_pipeline_tools_requirement.py

If the last output line is: [Takeda Pipeline Message] Go AHEAD! All the required tools are on the system., you're all set! If you see a tool missing, please install it using conda. A quick google search will give you the exact command for installation of any bioinformatics tool.

2.4 Create the Data Repository.
Please go to the directory variant-calling-pipeline using: cd variant-calling-pipeline ls

You'll see multiple folders such as: basics-stats, pipeline, data, simulate-data etc.

Note: You will not see the data directory the first time you're doing this process. Please make sure you're inside the variant-calling-pipeline directory and run the command: mkdir data

Please store your reference files and fastq files in this data repository.

Note: We usually do not store data on Bitbucket/GIT, and so you'll have to transfer the data using other means. We're still working on managing optimal sharing of data.

2.4 Naming of Data Files
The pipeline expects the data to be named in a certain way. Automating this process was a bit risky, because renaming a file programatically using python puts a risk of overwriting the file.

Step 1: Please go to the variant-calling-pipeline directory, and enter data directory:

cd variant-calling-pipeline/data
Step 2: Please create a directory for your cell lines, and enter there. You can use any name you want, specific to your run or data. Something like: cell-lines-amplicon-july-2020 would be more appropriate.

mkdir cell-lines-amplicon
cd cell-lines-amplicon
Step 3: Now you need to make different directories for all the cell lines you have:

mkdir HD-N62-F23 
mkdir HD-N62-F63
mkdir ZN23-F1
and so on...

Step 4: Please transfer the appropriate files associated into the Cell Line folders.

mv HD-N62-F23-file-1.fastq variant-calling-pipeline/data/cell-lines-amplicon/HD-N62-F23/HD-N62-F23_1.fastq
mv HD-N62-F23-file-2.fastq variant-calling-pipeline/data/cell-lines-amplicon/HD-N62-F23/HD-N62-F23_2.fastq
Since HD-N62-F23 (or any other cell line) will have two files (forward and reverse reads), please transfer both the files in the same cell line directory as shown above. Repeat this for all the cell lines: HD-N62-F63, ZN23-F1, ... etc

Note: It doesn't matter what your file names were, when you got them from your sequencing experiments. You have to make sure that the forward fastq file name is names as: cell_line_name_1.fastq and reverse fastq file is named as: cell_line_name_2.fastq.
Also, please make sure that the name: cell_line_name_1.fastq matches the directory name in which it exists.

For example: the files: HD-N62-F23_1.fastq and HD-N62-F23_2.fastq should be in direcotry HD-N62-F23

3. Running the Variant Calling Pipeline.
3.1 Standard Running & Testing the Pipeline
There are several aspects and ways in which you can run the pipeline. But the easiest way is using the following command:

python3 pipeline.py -r ../data/reference/processed/native.fasta -f1 ../data/fastq-files/_sample-name_/sample-name_1.fastq -f2 ../data/fastq-file/_sample-name_/sample-name_2.fastq

This will show a lot of text/messages on your screen and can take a lot of time to run depending upon the system you have.

3.2 Using the pipeline_runner.py Script
3.3 Running the application on the HPC cluster using PBS script.
A Portable Bash Script (PBS) is what you need to access Takeda's HPC. A script has already been written for you. The script can be submitted using the command: qsub pbs-script.pbs
