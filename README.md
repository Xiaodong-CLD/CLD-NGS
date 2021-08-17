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
You enter ```CUI3008@10.245.45.40``` in *Host Name* form field and click on *Open*. Enter you password and you should see something like:
		```cui3008@insight-hpc-login2  $```
        
Although, this pipeline can run on almost all the Unix/Linux environment, it is highly recommened to set up a *conda* environment on your Linux system.
In order to get *conda* on your system, go to [Miniconda's website](https://docs.conda.io/en/latest/miniconda.html). Right click on *Miniconda3 Linux 64-bit* for *Python 3* in the *Linux Installers* section; and select *Copy Link Address*.

Go back to you Putty window and write:
	```wget [CTRL+V  or  Right Click to Paste the Link]```
    
You should have successfully downloaded *conda* on your system. Now we will try to install *conda* on you HPC cluster account space. Please write:
	```bash Miniconda-Linux-64...```
  Please follow the instructions of the script and you will have to put in *Yes* in a couple of places.
  
Once you have completed this step, please close the Putty window and open Putty again and log into HPC again.

