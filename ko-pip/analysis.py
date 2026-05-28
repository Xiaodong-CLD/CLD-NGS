import argparse
import subprocess
import os
from pybedtools import BedTool

def run_command(command):
    """Execute the command (a string) using subprocess."""
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print(f"Error executing command: {command}")
        print(stderr.decode())
    return stdout.decode(), stderr.decode()

def align_contigs_to_reference(contigs_file, reference_file, output_prefix):
    """Align the contigs to the reference genome and return the path to the BAM file."""
    # Index the reference genome
    run_command(f"bwa index {reference_file}")

    # Align the contigs to the reference genome
    sam_file = f"{output_prefix}.sam"
    run_command(f"bwa mem {reference_file} {contigs_file} > {sam_file}")

    # Convert the SAM file to BAM
    bam_file = f"{output_prefix}.bam"
    run_command(f"samtools view -bS {sam_file} | samtools sort - > {bam_file}")

    # Index the BAM file
    run_command(f"samtools index {bam_file}")

    return bam_file

def compare_alignments(bam_files, output_file):
    """Compare the alignments from multiple samples."""
    # Convert the BAM files to BED
    bed_files = [f"{bam_file}.bed" for bam_file in bam_files]
    for bam_file, bed_file in zip(bam_files, bed_files):
        run_command(f"bedtools bamtobed -i {bam_file} > {bed_file}")

    # Use pybedtools to perform the intersection
    bed_tools = [BedTool(bed_file) for bed_file in bed_files]
    intersection = bed_tools[0].multi_intersect(i=bed_files[1:])
    intersection.saveas(output_file)

def find_uncovered_regions(reference_genome_file, alignment_bed_file, output_bed_file):
    """Find regions of the reference genome not covered by any contig."""
    # Create a BedTool for the reference genome
    genome = BedTool(reference_genome_file)

    # Create a BedTool for the alignment
    alignment = BedTool(alignment_bed_file)

    # Use the genome_coverage method to find regions of the genome that are not covered by any contig
    uncovered_regions = genome.genome_coverage(bg=True).subtract(alignment)

    # Save the uncovered regions to a file
    uncovered_regions.saveas(output_bed_file)

def perform_variant_analysis(contigs_file, reference_file, output_directory):
    """Perform variant analysis for SNPs, indels, and structural variants."""
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Your variant analysis commands go here
    # Example:
    variant_output = os.path.join(output_directory, "variant_analysis.txt")
    run_command(f"variant_analysis_command -i {contigs_file} -r {reference_file} -o {variant_output}")

def main(args):
    # Align the contigs to the reference genome
    bam_file = align_contigs_to_reference(args.contigs, args.reference, "alignment")

    # Compare the alignments
    compare_alignments([bam_file], os.path.join(args.output_dir, "alignment_comparison.bed"))

    # Find regions of the reference genome not covered by any contigs
    find_uncovered_regions(args.reference, os.path.join(args.output_dir, "alignment_comparison.bed"), os.path.join(args.output_dir, "uncovered_regions.bed"))

    # Perform variant analysis
    perform_variant_analysis(args.contigs, args.reference, args.output_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform analysis for a single sample.")
    parser.add_argument("-c", "--contigs", required=True, help="Path to the contigs FASTA file.")
    parser.add_argument("-r", "--reference", required=True, help="Path to the reference genome FASTA file.")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to the output directory.")
    args = parser.parse_args()

    main(args)
