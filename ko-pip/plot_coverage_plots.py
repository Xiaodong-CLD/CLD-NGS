#!/usr/bin/env python

'''
This script is a modification of the script located at:
	https://github.com/oneTakeda/voof/blob/master/amp-pip/pipeline.py
The above mentioned script on github is written by Tongyu Liu.
'''

import argparse
import os

def create_plotly_output_directory(input_dir):
    '''
    Create a directory called coverage-plot-output to store the graph for coverage.
    '''
    plotly_output_dir_path = os.path.join(input_dir, "coverage-plot-output")
    if os.path.exists(plotly_output_dir_path):
        return plotly_output_dir_path
    else:
        os.mkdir(plotly_output_dir_path)
        return plotly_output_dir_path


def read_consensus_file(consensus_file_path):
    '''
    We will read the consesnsus file produced by varscan to create a coverage plot with variants.
    '''
    with open(consensus_file_path) as f:
        raw = f.read()
    rows = [i.split('\t') for i in raw.split('\n') if i != '']
    consensus_rows = [[int(i[1]), int(i[14]) + int(i[15]), int(i[16]) + int(i[17])] for i in rows[1:]]

    return consensus_rows


# def plot_graphs(consensus_rows, cell_line_id, plotly_output_dir_path):
#     '''
#     For our graphs we will have two y-axes.
#     The first y-axis has the number of reads that supports the nucleotide suggested by the reference.
#     The second y-axis has the number of reads that do not support the nucleotide suggested by the reference, which means these reads support a variant nucleotide at this position.
#     '''
#     x_axis = [i[0] for i in consensus_rows]
#     y_axis_1 = [i[1] for i in consensus_rows]
#     y_axis_2 = [i[2] for i in consensus_rows]

#     try:
#         import plotly
#         import plotly.graph_objects as go
#     except ImportError:
#         print("[Takeda Pipeline Message] Plotly [A python package] is not installed. Please install it using: conda install -c plotly plotly ")
#         return

#     #Plot the graph.
#     fig = go.Figure()
#     fig.add_trace(go.Scatter(x = x_axis, y = y_axis_1, mode = 'lines', name = 'Reference'))
#     fig.add_trace(go.Scatter(x = x_axis, y = y_axis_2, mode = 'lines', name = 'Variant'))
#     fig.update_layout(title="Coverage for Reference & Variant Nucleotides for " + cell_line_id, xaxis_title = "Position", yaxis_title = "Number of Reads")

#     #Save the image.
#     image_path = os.path.join(plotly_output_dir_path, "coverage_plot.jpeg")
#     fig.write_image(image_path)

#     return True

def plot_graphs(consensus_rows, cell_line_id, plotly_output_dir_path):
    import matplotlib
    matplotlib.use("Agg")               # no display needed
    import matplotlib.pyplot as plt

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
    fig.savefig(os.path.join(plotly_output_dir_path, "coverage_plot.png"), dpi=150)
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
        print("[Takeda Pipeline Message] Sorry, the coverage file is not present. Please make sure a file like: varscan-output/HD-N7-F12.consensus.vcft")
        return False

    input_dir = '/'.join(input_dir.split('/')[:-1])
    input_dir = os.path.dirname(input_dir)

    #Create a output direcotry for plotly.
    plotly_output_dir_path = create_plotly_output_directory(input_dir)

    consensus_rows = read_consensus_file(consensus_file_path)
    plot_graphs(consensus_rows, cell_line_id, plotly_output_dir_path)


if __name__ == "__main__":
    plot_coverage_plots_manager(None)
