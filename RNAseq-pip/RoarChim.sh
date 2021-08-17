#!/bin/bash

set -e

# source /mnt/efs/home/oqh2473/miniconda3/bin/activate
# conda activate RNAseq

while getopts d:g:o: flag
do
    case "${flag}" in
        d) fq_dir=${OPTARG};;
        g) star_genome=${OPTARG};;
        o) out_dir=${OPTARG};;
    esac
done

# These four varibles need to be edited to adapt to other gene of interest.
# Please also scan this script to replace "ADAMTS13" to gene of interest in future usage.
gtf=/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/ADAMTS13.gtf
ref_seq=/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/adam_cho.fa
goi_seq=/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/ADAMTS13.fa
bed=/mnt/efs/home/cui3008/reference_genome/RNAseq_genome/ncbi/ADAMTS13.bed


base=`basename $fq_dir`
fq_files=($fq_dir/*)
fq_file="${fq_files[0]}"
extension="${fq_file##*.}"

echo "Takeda pipeline: Processing Sample $base"
if [ ! -s $out_dir/$base ]
then
  mkdir $out_dir/$base
else echo "Takeda pipeline: output directory for this sample already exsit"
fi

out_dir=$out_dir/$base

# fastp to run quality control
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

# Run alignment and chimeric junction detection
cd $fq_dir
if [ ! -s $base.Chimeric.out.junction ]
then
  echo "Takeda pipeline: Starting STAR alignment and generating chimeric junctions..."
  STAR --runMode alignReads --twopassMode Basic --runThreadN 12 --genomeDir $star_genome \
  --readFilesIn $fq_dir/*1.trimmed* $fq_dir/*2.trimmed*  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $base. --sjdbOverhang 149 \
  --outSAMunmapped Within --chimSegmentMin 10 --chimOutType Junctions --chimFilter banGenomicN --chimMainSegmentMultNmax 1 --chimOutJunctionFormat 1 --chimMultimapNmax 1
  cat $base.Chimeric.out.junction |  awk '/ADAMTS13/{ print }' > $out_dir/$base.adam.chimeric.txt
else echo "Takeda pipeline: STAR already performed... skipping STAR alignment step"
fi

# Marker duplicates
cd $fq_dir
if [ ! -s $base.sorted.bam ]
then
  echo "Takeda pipeline: sorting bam..."
  samtools sort -O BAM -o $base.sorted.bam -T $base.bam.sort *.sortedByCoord.out.bam
else echo "Takeda pipeline: bam file sorted... skipping"
fi

# mark duplicates, need java 8 to execute#
cd $fq_dir
if [ ! -s $base.md.bam ]
then
  echo "Takeda pipeline: Marking duplicates..."
  picard MarkDuplicates I=$base.sorted.bam O=$base.md.bam ASSUME_SORTED=true METRICS_FILE=$base.md.bam.markdup.metrics VALIDATION_STRINGENCY=SILENT
else echo "Takeda pipeline: marked duplicates bam file already exits...skipping"
fi

# samtools index, filter and pileup
cd $fq_dir
if [ ! -s *.hqa.mpileup ]
then
  echo "Takeda pipeline: extracting high quality alignment and removing duplicates..."
  samtools index $base.md.bam
  samtools view -b -h -f 3 -F 4 -F 256 -F 1024 -F 2048 -q 30 $base.md.bam  > $base.hqa.bam
  samtools index $fq_dir/$base.hqa.bam
  samtools mpileup -f  $ref_seq  $base.hqa.bam -l $bed -A > $base.targeted.hqa.mpileup
else echo "Takeda pipeline: high quality alignment (hqa) bam already exits...skipping"
fi

# varscan to call variant
cd $fq_dir
if [ ! -s $out_dir/*snp.vcf ]
then
  varscan pileup2snp $base.targeted.hqa.mpileup --min-coverage 300 > $fq_dir/$base.targeted.hqa.mpileup.snp.vcf
  varscan pileup2indel $base.targeted.hqa.mpileup  > $fq_dir/$base.targeted.hqa.mpileup.indel.vcf
  varscan pileup2cns $base.targeted.hqa.mpileup  > $fq_dir/$base.targeted.hqa.mpileup.consensus.vcf
  mv  *vcf $out_dir/
else echo "Takeda pipeline: vcf files already exist...skipping"
fi

# pindel, I use Shivam's pindel function to run it.
cd $fq_dir
adam=$(cat $bed |  awk '{print $1 ":" $2+1 "-" $3}')
samtools view -b $base.hqa.bam $adam > $base.targeted.hqa.bam
samtools index $base.targeted.hqa.bam
if [ ! -s $out_dir/pindel-output/*Summary* ]
then
  python /mnt/efs/home/cui3008/tongyu/bin/Run_pindel_only.py \
  -r $goi_seq -b $base.targeted.hqa.bam
  mv $fq_dir/pindel-output $out_dir/$base.pindel
else echo "Takeda pipeline: pindel outputs already exist...skipping"
fi

# htseq count gene of interest only
if [ ! -s $out_dir/$base.htseq ]
then
  htseq-count -f bam -t exon -i gene_id --stranded=no $fq_dir/$base.Aligned.sortedByCoord.out.bam $gtf > $fq_dir/$base.htseq
  mv $fq_dir/$base.htseq $out_dir/$base.htseq
else echo "Takeda pipeline: htseq files already exist...skipping"
fi

# calculate chimeric reads percentage
# four arguments for chimeric_stat.R as input

if [ -s $out_dir/$base.adam.chimeric.txt -a -s $out_dir/$base.htseq  ]
then
  echo 'Takeda pipeline: Summarizing chimeric transcripts counts...'
  Rscript /mnt/efs/home/cui3008/tongyu/bin/chimeric_stat.R \
  $out_dir/$base.adam.chimeric.txt $out_dir/$base.htseq $base $out_dir
else echo 'Takeda pipeline: missing inputs for chimeric transcripts report.'
fi

# extract chimeric reads from fastq files
if [ -s $out_dir/$base.chimeric.reads.name.txt ]
then
  echo  'Takeda pipeline: Extracting chimeric reads...'
  /mnt/efs/home/cui3008/tongyu/bin/filterbyname.sh in=$fq_dir/*1.trimmed* in2=$fq_dir/*2.trimmed* \
  out=$out_dir/$base.fused.1.fq out2=$out_dir/$base.fused.1.fq names=$out_dir/$base.chimeric.reads.name.txt include=t
else echo 'Takeda pipeline: no chimeric reads to extract...skipping'
fi

# make coverage plot #
if [ -s $out_dir/$base.targeted.hqa.mpileup.consensus.vcf  ]
then
  echo 'Takeda pipeline: Generating coverage plots...'
  Rscript /mnt/efs/home/cui3008/tongyu/bin/coverage_plot.R \
  $out_dir/$base.targeted.hqa.mpileup.consensus.vcf $base $out_dir
else echo 'Takeda pipeline: missing inputs for coverage plot'
fi
# final information
echo 'Takeda pipeline: RoarChim pipeline finished. Please check final results in output folder'
