import subprocess
import os
import random

def check_tools():
    '''
    Before running the pipeline, we need to make sure if all the requirements are installed on the system.
    Please run this script to check if all the tools are installed.

    If you need to add a new tool to the list of tools required, please add that in the follwoing python list variable as a string.
    '''
    tools_list = ['samtools', 'bamtools', 'fastp', 'bwa', 'varscan', 'pindel', 'pindel2vcf']
    all_ok = True
    for tool in tools_list:
        try:
            bash_output = subprocess.call([tool], stdout = subprocess.PIPE, stderr = subprocess.PIPE) 
            print("[Takeda Pipeline Message] Success, system has {} installed on it".format(tool))

        except (subprocess.CalledProcessError, FileNotFoundError):
            print("[Takeda Pipeline Message] A Particular Bioinformatics tools is missing, please install it using conda.")
            print("[Takeda Pipeline Message] Tool Name: {}".format(tool))
            all_ok = False
            
    #See if plotly is installed on the system.
    try:
        import plotly
    except ImportError:
         print("[Takeda Pipeline Message] A Particular Bioinformatics tools is missing, please install it using conda.")
         print("[Takeda Pipeline Message] Tool Name: Plotly Python Module")

    if all_ok:
        print("[Takeda Pipeline Message] Go AHEAD! All the required tools are on the system.")
    else:
        print("[Takeda Pipeline Message] Please install the required tools before running the pipeline.")

        
if __name__ == "__main__":
	check_tools()
	pass
