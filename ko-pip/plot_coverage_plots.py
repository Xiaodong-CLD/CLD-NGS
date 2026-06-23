#!/usr/bin/env python

'''
This script is a modification of the script located at:
	https://github.com/oneTakeda/voof/blob/master/amp-pip/pipeline.py
The above mentioned script on github is written by Tongyu Liu.

The script was modified by Cindy Ren 2026-06-02. Replaced the plotly-based
`plot_graphs` implementation with a matplotlib version that writes a PNG 
instead of a JPEG, removing the plotly/kaleido dependency.
'''

import argparse
import os
import matplotlib
matplotlib.use("Agg") # no display needed
import matplotlib.pyplot as plt

def create_plots_output_directory(input_dir):
    '''
    Create a directory called plots to store the coverage plot.
    '''
    plots_output_dir_path = os.path.join(input_dir, "plots")
    if os.path.exists(plots_output_dir_path):
        return plots_output_dir_path
    else:
        os.mkdir(plots_output_dir_path)
        return plots_output_dir_path


def _safe_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def read_consensus_file(consensus_file_path):
    '''
    Read the VarScan `mpileup2cns` consensus table and return
    [position, ref_reads, var_reads] per position for the coverage plot.

    mpileup2cns uses the packed genotype schema: ref/var read counts come from
    the `Cons:Cov:Reads1:Reads2:Freq:P-value` field at col 4 (Reads1 / Reads2);
    `-` placeholders at uncovered positions are read as 0.
    '''
    with open(consensus_file_path) as f:
        raw = f.read()
    rows = [i.split('\t') for i in raw.split('\n') if i != '']
    if not rows:
        return []
    consensus_rows = []
    for i in rows[1:]:
        packed = i[4].split(':') if len(i) > 4 else []
        reads1 = _safe_int(packed[2]) if len(packed) > 2 else 0
        reads2 = _safe_int(packed[3]) if len(packed) > 3 else 0
        consensus_rows.append([int(i[1]), reads1, reads2])

    return consensus_rows


def plot_graphs(consensus_rows, cell_line_id, plots_output_dir_path):
    x = [i[0] for i in consensus_rows]
    y_ref = [i[1] for i in consensus_rows]
    y_var = [i[2] for i in consensus_rows]

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(x, y_ref, label="Reference", linewidth=1)
    ax.plot(x, y_var, label="Variant",   linewidth=1, color="crimson")
    ax.set_xlabel("Position")
    ax.set_ylabel("Number of Reads")
    ax.set_title(f"Coverage for Reference & Variant Nucleotides — {cell_line_id}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(plots_output_dir_path, "coverage_plot.png"), dpi=150)
    plt.close(fig)
    return True


def plot_coverage_plots_manager(input_dir):
    if input_dir is None:
        #Take in the reference file as arguments.
        parser = argparse.ArgumentParser()

        #Add argument options
        parser.add_argument("-i", "--input_dir", help="Path of directory which has fastq files (paired).", required = True)

        #Fetch the argument data.
        args = vars(parser.parse_args())
        input_dir = args['input_dir']

    #Get the consensus file.
    cell_line_id = input_dir.rstrip('/').split('/')[-1]
    # cell_line_id:Spike-Mut10-1_S7_L001.consensus.vcf
    print("[Takeda Pipeline Message] Creating coverage plot for cell line: {} using plotly".format(cell_line_id))

    consensus_file_path = input_dir
    if not os.path.exists(consensus_file_path):
        print("[Takeda Pipeline Message] Sorry, the coverage file is not present. Please make sure a file like: snv/HD-N7-F12.consensus.vcf is present.")
        return False

    input_dir = '/'.join(input_dir.split('/')[:-1])
    input_dir = os.path.dirname(input_dir)

    #Create a output direcotry for the coverage plot.
    plots_output_dir_path = create_plots_output_directory(input_dir)

    consensus_rows = read_consensus_file(consensus_file_path)
    plot_graphs(consensus_rows, cell_line_id, plots_output_dir_path)


if __name__ == "__main__":
    plot_coverage_plots_manager(None)
