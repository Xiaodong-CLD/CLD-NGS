
# Cell Line Development amplicon-seq pipeline

## User Manual Amp Var Calling.pdf contains the details about running this pipeline, including setting up environment, preparing data, and modifying code to fit in personal environment. Below is the short version of instruciton on running this pipeline. 

## Setting up working environment
1. Setup conda installation following the instructions in the main repository README or follow the user manual for this pipeline.
2. `conda activate variant-calling`. `variant-calling.yml` contains an export of the packages in case conda needs to be reinstalled.

## Running pipeline
### Preparing configuration spreadsheet
This version of the amplicon pipeline requires a configuration spreadsheet to run. Please see `example.xlsx` for an example spreadsheet. The speadsheet keeps track of each sample (one per line), and the reference sequence and fastq filenames for each. The `Clone ID` column is optional, and matches replicates together under a Clone ID.

### Directory structure
All fastq files go in a `fastq/` directory, and all GOI `.fasta` files go in a `reference/` directory. The configuration spreadsheet goes in the root directory.

	<input directory>/
	├─ Fastq/
	│  ├─ All forward and reverse fastq files referenced in spreadsheet go here
	├─ reference/
	│  ├─ All reference .fasta files in spreadsheet go here
	├─ <run name>.xlsx
	├─ job.script

### Usage and options
	usage: pipeline_combined.py [-h] input_dir

	positional arguments:
		input_dir   Path to input directory

### Starting the run
See `job.script` for an example Slurm script.
Submit the job to the queue with `sbatch job.script`.

### Outputs
- A folder with each sample ID will be created that contains the outputs for that sample.
- All varscan and pindel reports are contained in each sample folder. For the varscan snp report, the amino acid annotation is also included as an additional file.
- A summary varscan snp report (named `<run name>.snp.txt`) is written to the root directory of the run. If the `Clone ID` column in the configuration spreadsheet is filled for all samples, this summary snp report will consider the replicates. If any sample is missing the `Clone ID` field, the summary snp report is split out by sample.
