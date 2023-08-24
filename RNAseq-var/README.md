# Cell Line Development RNAseq-var pipeline
### User Manual Rna Var Calling.pdf contains the details about running this pipeline, including setting up environment, preparing data, and modifying code to fit in personal environment. Below is the short version of instruciton on running this pipeline. 


### Setting up the working environment.
1. Setup conda installation following the instructions in the main repository README.
2. `conda activate RNAseq`. `RNAseq.yml` contains an export of the packages in case conda needs to be reinstalled.

### Directory structure
All fastq files go in a `Fastq/` directory, all GOI `.fasta` files go in a `goi/` directory, the reference genome `.fasta` file goes into `reference/`. The configuration spreadsheet goes in the root directory.

	<input directory>/
	├─ Fastq/
	│  ├─ All forward and reverse fastq files referenced in spreadsheet go here
	├─ reference/
	│  ├─ Reference .fasta file in spreadsheet goes here
	├─ goi/
	│  ├─ All goi .fasta files in spreadsheet goes here
	├─ <run name>.xlsx
	├─ star_job.script
	├─ job.script

### Building Star index
Modify the file path in the `star_job.script` and submit a job to the queue with `sbatch star_job.script`. Wait until this job is finished and continue to perform RNAseq variant calling. 
#### Outputs
This `star_job.script` will automatically concatenate goi and reference sequences. And build star index for goi, reference, and combined sequences. Once it is completed, in the dirctory there would be 3 star index directories:ref_star_dir, goi_star_dir, combined_star_dir.

### Running the Variant Calling Pipeline.
See `job.script` for an example Slurm script.
Submit the job to the queue with `sbatch job.script`.

#### Outputs
- A folder with each sample ID will be created that contains the outputs for that sample.
- All varscan and pindel reports are contained in each sample folder. For the varscan snp report, the amino acid annotation is also included as an additional file.
- A summary varscan snp report (named `<run name>.snp.txt`) is written to the root directory of the run. If the `Clone ID` column in the configuration spreadsheet is filled for all samples, this summary snp report will consider the replicates. If any sample is missing the `Clone ID` field, the summary snp report is split out by sample.
	


	
