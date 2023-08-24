#!/bin/bash
set -e

while getopts d:o:c:r:g:s: flag
do
    case "${flag}" in
        d) fq_dir=${OPTARG};;
        o) out_dir=${OPTARG};;
        c) combined_dir=${OPTARG};;
        r) ref_dir=${OPTARG};;
        g) goi_dir=${OPTARG};;
        s) goi_seq=${OPTARG};;
    esac
done

base=`basename $fq_dir`
fq_files=($fq_dir/*)
fq_file="${fq_files[0]}"
extension="${fq_file##*.}"

echo "[Takeda Pipeline Message] Processing Sample $base"
if [ ! -s $out_dir/$base ]
then
  mkdir $out_dir/$base
else echo "[Takeda Pipeline Message] output directory for this sample already exists"
fi

out_dir=$out_dir/$base

# Fastp to run quality control
cd $fq_dir
if [ ! -s *1.trimmed* ]
then
  if [ "$extension" = "fastq" ]
  then
    echo '[Takeda Pipeline Message] Input is fastq, running fastp...'
    fastp --in1  *1.fastq  --in2  *2.fastq  --out1 $fq_dir/$base.R1.trimmed.fastq  --out2 $fq_dir/$base.R2.trimmed.fastq  \
    -w 16 -j  $fq_dir/$base.fastp.json  -h  $fq_dir/$base.fastp.html
  elif [ "$extension" = "gz" ]
  then
    echo '[Takeda Pipeline Message] Input is .fastq.gz, running fastp...'
    fastp --in1  *R1*.gz  --in2  *R2*.gz  --out1 $fq_dir/$base.R1.trimmed.fq.gz  --out2 $fq_dir/$base.R2.trimmed.fq.gz  \
    -w 16 -j  $fq_dir/$base.fastp.json  -h  $fq_dir/$base.fastp.html
  else echo '[Takeda Pipeline Message] Please check input format, either fastq or .fq.gz'
  fi
else echo "[Takeda Pipeline Message] Trimmed fastq already exist... skipping QC step"
fi

# Count number of reads
echo "[Takeda Pipeline Message] Counting the number of reads after fastp..."
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

# Run alignment to human reference
if [ ! -s $base.human.bam ]
then
  echo "[Takeda Pipeline Message] Starting STAR alignment to human reference..."
  STAR --genomeDir $ref_dir --readFilesIn $fq_dir/*1.trimmed* $fq_dir/*2.trimmed* --readFilesCommand zcat --runThreadN 24 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $base.human.
  # Rename output
  mv $base.human.Aligned.sortedByCoord.out.bam $base.human.bam

  echo "[Takeda Pipeline Massage] Counting the reads mapped to human"
  # Count reads mapped to human
  mapped_human_reads=$(samtools view -c -F 4 $base.human.bam)
  echo "Reads mapped to human reference when aligning separately: $mapped_human_reads" >> $out_dir/read_counts.txt
fi

# Run alignment to goi reference
# To sparate from below goi from combined analysis, use Goi here
if [ ! -s $base.Goi.bam ]
then
  echo "[Takeda Pipeline Message] Starting STAR alignment to goi reference..."
  STAR --genomeDir $goi_dir --readFilesIn $fq_dir/*1.trimmed* $fq_dir/*2.trimmed* --readFilesCommand zcat --runThreadN 24 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $base.Goi.
  # Rename output
  mv $base.Goi.Aligned.sortedByCoord.out.bam $base.Goi.bam

  echo "[Takeda Pipeline Massage] Counting the reads mapped to goi"
  # Count reads mapped to goi
  mapped_goi_reads=$(samtools view -c -F 4 $base.Goi.bam)
  echo "Reads mapped to goi when aligning separately: $mapped_goi_reads" >> $out_dir/read_counts.txt
fi


# Run alignment with STAR to combined reference
if [ ! -s $base.combined.Aligned.out.bam ]
then
  echo "[Takeda Pipeline Message] Starting STAR alignment to combined reference..."
  STAR --genomeDir $combined_dir --readFilesIn $fq_dir/*1.trimmed* $fq_dir/*2.trimmed* --readFilesCommand zcat --runThreadN 24 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $base.combined.
  
  # Rename
  mv $base.combined.Aligned.sortedByCoord.out.bam $base.combined.Aligned.out.bam

  # Index BAM  file
  echo "[Takeda Pipeline Message] Index BAM file..."
  samtools index $base.combined.Aligned.out.bam
  
  # Filter out reads mapped to human reference and keep only reads mapped to GOI
  GOI_reference_name=$(grep -m1 ">" $goi_seq | cut -d' ' -f1 | tr -d ">")
  samtools view -b $base.combined.Aligned.out.bam "$GOI_reference_name" > $base.goi.bam
  
  # Count total reads, reads mapped to goi, unmapped reads, and reads mapped to human
  total_reads=$(samtools view -c $base.combined.Aligned.out.bam)
  unmapped_reads=$(samtools view -c -f 4 $base.combined.Aligned.out.bam)
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
  echo "[Takeda Pipeline Message] Sort and index GOI bam file..."
  samtools sort $base.goi.bam -o $base.goi.sorted.bam
  samtools index $base.goi.sorted.bam
else echo "[Takeda Pipeline Message] GOI sorted bam file already exists... skipping"
fi

# Create mpileup file
cd $fq_dir
if [ ! -s $base.goi.mpileup ]
then
  echo "[Takeda Pipeline Message] Creating mpileup file..."
  samtools mpileup -f $goi_seq $base.goi.sorted.bam -o $base.goi.mpileup 
else 
  echo "[Takeda Pipeline Message] Mpileup file already exists... skipping"
fi

############################################################
# Variant calling using VarScan
cd $fq_dir
if [ ! -d $out_dir/$base.snp ]
then
  echo "[Takeda Pipeline Message] Starting VarScan..."
  varscan pileup2snp $base.goi.mpileup --min-coverage 300 > $fq_dir/$base.goi.mpileup.snp.vcf
  varscan pileup2indel $base.goi.mpileup  > $fq_dir/$base.goi.mpileup.indel.vcf
  varscan pileup2cns $base.goi.mpileup  > $fq_dir/$base.goi.mpileup.consensus.vcf

  # Move VCF files to the output directory
  echo "[Takeda Pipeline Message] Moving vcf files to output directory..."
  if [ ! -s $out_dir/$base.snp ]
  then
    mkdir $out_dir/$base.snp
  else echo "[Takeda Pipeline Message] output directory for snp files"
  fi
  mv  *vcf $out_dir/$base.snp
else echo "[Takeda Pipeline Message] VCF files already exist... skipping VarScan"
fi

##############################################################
# Pindel
cd $fq_dir
if [ ! -d $out_dir/$base.pindel ]
then
  echo "[Takeda Pipeline Message] Running Pindel..."
  python /mnt/efs/home/zds0975/pipeline_test/voof/RNAseq-var/bin/Run_pindel_only.py -r $goi_seq -b $base.goi.sorted.bam
  
  # Move pindel files to the output directory
  echo "[Takeda Pipeline Message] Moving pindel files to output directory..."
  mv $fq_dir/pindel-output $out_dir/$base.pindel
else 
  echo "[Takeda Pipeline Message] Pindel outputs already exist... skipping"
fi


# Coverage plot
cd $fq_dir
cd $out_dir/$base.snp
if [ -s $base.goi.mpileup.consensus.vcf ]
then
  echo '[Takeda Pipeline Message] Generating coverage plots...'
  Rscript /mnt/efs/home/zds0975/pipeline_test/voof/RNAseq-var/bin/coverage_plot.R \
  $out_dir/$base.snp/$base.goi.mpileup.consensus.vcf $base $out_dir
else echo '[Takeda Pipeline Message] missing inputs for coverage plot'
fi

# final information
echo '[Takeda Pipeline Message] RNAseq variant-calling pipeline finished. Please check final results in output folder'
