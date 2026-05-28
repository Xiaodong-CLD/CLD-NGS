#!/bin/bash

set -e

while getopts d:o:c:r:g: flag
do
    case "${flag}" in
        d) fq_dir=${OPTARG};;
        o) out_dir=${OPTARG};;
        c) combined_seq=${OPTARG};;
        r) ref_seq=${OPTARG};;
        g) goi_seq=${OPTARG};;
    esac
done

base=`basename $fq_dir`
fq_files=($fq_dir/*)
fq_file="${fq_files[0]}"
extension="${fq_file##*.}"

echo "Takeda pipeline: Processing Sample $base"
if [ ! -s $out_dir/$base ]
then
  mkdir $out_dir/$base
else echo "Takeda pipeline: output directory for this sample already exists"
fi

out_dir=$out_dir/$base

# Fastp to run quality control
cd $fq_dir
if [ ! -s *1.trimmed* ]
then
  if [ "$extension" = "fastq" ]
  then
    echo 'Takeda pipeline: Input is fastq, running fastp...'
    fastp --in1  *1.fastq  --in2  *2.fastq  --out1 $fq_dir/$base.R1.trimmed.fastq  --out2 $fq_dir/$base.R2.trimmed.fastq  \
    -w 16 -j  $fq_dir/$base.fastp.json  -h  $fq_dir/$base.fastp.html
  elif [ "$extension" = "gz" ]
  then
    echo 'Takeda pipeline: Input is .fastq.gz, running fastp...'
    fastp --in1  *R1*.gz  --in2  *R2*.gz  --out1 $fq_dir/$base.R1.trimmed.fq.gz  --out2 $fq_dir/$base.R2.trimmed.fq.gz  \
    -w 16 -j  $fq_dir/$base.fastp.json  -h  $fq_dir/$base.fastp.html
  else echo 'Takeda pipeline: Please check input format, either fastq or .fq.gz'
  fi
else echo "Takeda pipeline: Trimmed fastq already exist... skipping QC step"
fi

# Count number of reads
echo "Takeda pipeline: Counting the number of reads after fastp..."
if [ "$extension" = "fastq" ]
then
  count_fastp_R1=$(grep -c "^+$" $fq_dir/$base.R1.trimmed.fastq)
  count_fastp_R2=$(grep -c "^+$" $fq_dir/$base.R2.trimmed.fastq)
elif [ "$extension" = "gz" ]
then
  count_fastp_R1=$(zgrep -c "^+$" $fq_dir/$base.R1.trimmed.fq.gz)
  count_fastp_R2=$(zgrep -c "^+$" $fq_dir/$base.R2.trimmed.fq.gz)
fi
total_fastp_reads=$((count_fastp_R1 + count_fastp_R2))
echo "Total reads after fastp: $total_fastp_reads" >> $out_dir/read_counts.txt

###################################
# cd $fq_dir
# # Create a combined reference
# if [ ! -s combined.fasta ]
# then
#   echo "Takeda pipeline: Creating a combined reference..."
#   cat $ref_seq $goi_seq > combined.fasta
#   bwa index combined.fasta
# fi
##########################################################
# Run alignment to human reference
if [ ! -s $base.human.bam ]
then
  echo "Takeda pipeline: Starting BWA alignment to human reference..."
  bwa mem $ref_seq $fq_dir/*1.trimmed* $fq_dir/*2.trimmed* | samtools view -S -b > $base.human.bam

  # Count reads mapped to human
  mapped_human_reads=$(samtools view -c -F 4 $base.human.bam)
  echo "Reads mapped to human reference when aligning separately: $mapped_human_reads" >> $out_dir/read_counts.txt
fi

# Run alignment to goi reference
if [ ! -s $base.Goi.bam ]
then
  echo "Takeda pipeline: Starting BWA alignment to goi reference..."
  bwa mem $goi_seq $fq_dir/*1.trimmed* $fq_dir/*2.trimmed* | samtools view -S -b > $base.Goi.bam

  # Count reads mapped to goi
  mapped_goi_reads=$(samtools view -c -F 4 $base.Goi.bam)
  echo "Reads mapped to goi when aligning separately: $mapped_goi_reads" >> $out_dir/read_counts.txt
fi


# Run alignment to the combined reference
if [ ! -s $base.combined.bam ]
then
  echo "Takeda pipeline: Starting BWA alignment to combined reference..."
  bwa mem $combined_seq $fq_dir/*1.trimmed* $fq_dir/*2.trimmed* | samtools view -S -b > $base.combined.bam
  
  # Sort and index the BAM file
  samtools sort $base.combined.bam -o $base.combined.sorted.bam 
  samtools index $base.combined.sorted.bam

  # Filter out reads mapped to human reference and keep only reads mapped to GOI
  GOI_reference_name=$(grep -m1 ">" $goi_seq | cut -d' ' -f1 | tr -d ">")
  samtools view -b $base.combined.sorted.bam "$GOI_reference_name" > $base.goi.bam
  
  # Count total reads, reads mapped to goi, unmapped reads, and reads mapped to human
  total_reads=$(samtools view -c $base.combined.bam)
  unmapped_reads=$(samtools view -c -f 4 $base.combined.bam)
  mapped_to_goi_reads=$(samtools view -c -F 4 $base.goi.bam)
  mapped_to_human_reads=$((total_reads - mapped_to_goi_reads - unmapped_reads))

  echo "Total reads after alignment to combined reference: $total_reads" >> $out_dir/read_counts.txt
  echo "Reads mapped to goi when aligning simultaneously: $mapped_to_goi_reads" >> $out_dir/read_counts.txt
  echo "Reads mapped to human when aligning simultaneously: $mapped_to_human_reads" >> $out_dir/read_counts.txt
fi


# Sort and index goi bam file
cd $fq_dir
if [ ! -s $base.goi.sorted.bam ]
then
  echo "Takeda pipeline: Sort and index GOI bam file..."
  samtools sort $base.goi.bam -o $base.goi.sorted.bam
  samtools index $base.goi.sorted.bam
else echo "Takeda pipeline: GOI sorted bam file already exists... skipping"
fi

# Create mpileup file
cd $fq_dir
if [ ! -s $base.goi.mpileup ]
then
  echo "Takeda pipeline: Creating mpileup file..."
  samtools mpileup -f $goi_seq $base.goi.sorted.bam -o $base.goi.mpileup 
else 
  echo "Takeda pipeline: Mpileup file already exists... skipping"
fi

############################################################
# Variant calling using VarScan
cd $fq_dir
if [ ! -d $out_dir/$base.snp ]
then
  echo "Takeda pipeline: Starting VarScan..."
  varscan pileup2snp $base.goi.mpileup --min-coverage 300 > $fq_dir/$base.goi.mpileup.snp.vcf
  varscan pileup2indel $base.goi.mpileup  > $fq_dir/$base.goi.mpileup.indel.vcf
  varscan pileup2cns $base.goi.mpileup  > $fq_dir/$base.goi.mpileup.consensus.vcf

  # Amino acid annotation on VarScan snp report
  # echo "Takeda pipeline: Annotating variants..."
  # python3 /mnt/efs/home/zds0975/pipeline_test/voof/amino-anno/from_varscan.py  $fq_dir/$base.goi.mpileup.snp.vcf $goi_seq -o $fq_dir

  # Move VCF files to the output directory
  echo "Takeda pipeline: Moving vcf files to output directory..."
  if [ ! -s $out_dir/$base.snp ]
  then
    mkdir $out_dir/$base.snp
  else echo "Takeda pipeline: output directory for snp files"
  fi
  mv  *vcf $out_dir/$base.snp
else echo "Takeda pipeline: VCF files already exist... skipping VarScan"
fi

##############################################################
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# Pindel
cd $fq_dir
if [ ! -d $out_dir/$base.pindel ]
then
  echo "Takeda pipeline: Running Pindel..."
  python "$SCRIPT_DIR/bin/Run_pindel_only.py" -r $goi_seq -b $base.goi.sorted.bam
  
  # Move pindel files to the output directory
  echo "Takeda pipeline: Moving pindel files to output directory..."
  mv $fq_dir/pindel-output $out_dir/$base.pindel
else 
  echo "Takeda pipeline: Pindel outputs already exist... skipping"
fi


# Coverage plot
cd $fq_dir
cd $out_dir/$base.snp
if [ -s $base.goi.mpileup.consensus.vcf ]
then
  echo 'Takeda pipeline: Generating coverage plots...'
  Rscript "$SCRIPT_DIR/bin/coverage_plot.R" \
  $out_dir/$base.snp/$base.goi.mpileup.consensus.vcf $base $out_dir
else echo 'Takeda pipeline: missing inputs for coverage plot'
fi

# final information
echo 'Takeda pipeline: RNAseq variant-calling pipeline finished. Please check final results in output folder'
