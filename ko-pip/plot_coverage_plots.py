#!/usr/bin/env python

'''
This script is a modification of the script located at:
	https://github.com/oneTakeda/voof/blob/master/amp-pip/pipeline.py
The above mentioned script on github is written by Tongyu Liu.

Coverage plot for the KO pipeline. Reads the VarScan ``mpileup2cns`` consensus
table and draws reference- vs variant-supporting read depth across the whole
amplicon.

Pinned to plotly 5.x and kaleido 0.2.1 to avoid the plotly 6.x / kaleido 0.3.x
breakage (chromium not bundled, needs system Chrome).
'''

import argparse
import os


def create_plots_output_directory(input_dir):
    '''
    Create a directory called plots to store the coverage plot.
    '''
    plots_output_dir_path = os.path.join(input_dir, "plots")
    os.makedirs(plots_output_dir_path, exist_ok=True)
    return plots_output_dir_path


def _safe_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def get_reference_length(reference_file_path):
    '''
    Return the total sequence length (bp) of the reference FASTA. Used to pad
    the coverage with zeros out to the end of the amplicon after the last
    covered position.
    '''
    length = 0
    with open(reference_file_path) as f:
        for line in f:
            if not line.startswith('>'):
                length += len(line.strip())
    return length


def read_consensus_file(consensus_file_path, reference_length=None):
    '''
    Read the VarScan ``mpileup2cns`` consensus table and return
    ``[position, ref_reads, var_reads]`` per position for the coverage plot.

    The consensus file is the plain tabular ``mpileup2cns`` output (despite the
    ``.consensus.vcf`` name): column 4 (``Reads1``) is the reference-supporting
    read count and column 5 (``Reads2``) is the variant-supporting read count.

    Positions with no consensus entry are padded with zero coverage so the plot
    spans the full amplicon instead of starting partway in: positions before the
    first covered position, and -- when ``reference_length`` is given -- the
    positions after the last covered one up to the reference length.
    '''
    with open(consensus_file_path) as f:
        raw = f.read()
    rows = [i.split('\t') for i in raw.split('\n') if i != '']
    if not rows:
        return []

    consensus_rows = []
    for i in rows[1:]:
        if len(i) < 6:
            continue
        consensus_rows.append([_safe_int(i[1]), _safe_int(i[4]), _safe_int(i[5])])

    if not consensus_rows:
        return consensus_rows

    consensus_rows.sort(key=lambda r: r[0])
    first_pos = consensus_rows[0][0]
    last_pos = consensus_rows[-1][0]
    end = reference_length if reference_length is not None else last_pos

    padded = [[pos, 0, 0] for pos in range(1, first_pos)]
    padded += consensus_rows
    padded += [[pos, 0, 0] for pos in range(last_pos + 1, end + 1)]

    return padded


def plot_graphs(consensus_rows, cell_line_id, plots_output_dir_path):
    '''
    Two traces per plot:
      - Reference: reads supporting the reference nucleotide at each position.
      - Variant:   reads supporting a non-reference (variant) nucleotide.

    Writes a static ``coverage_plot.jpeg`` via kaleido (pinned to 0.2.1, which
    bundles chromium and needs no system Chrome).
    '''
    x_axis = [i[0] for i in consensus_rows]
    y_axis_ref = [i[1] for i in consensus_rows]
    y_axis_var = [i[2] for i in consensus_rows]

    try:
        import plotly.graph_objects as go
    except ImportError:
        print("[Takeda Pipeline Message] Plotly [A python package] is not installed. "
              "Please install it using: conda install -c plotly plotly")
        return False

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x_axis, y=y_axis_ref, mode='lines', name='Reference'))
    fig.add_trace(go.Scatter(x=x_axis, y=y_axis_var, mode='lines', name='Variant'))
    fig.update_layout(
        title="Coverage for Reference & Variant Nucleotides for " + cell_line_id,
        xaxis_title="Position",
        yaxis_title="Number of Reads",
    )

    image_path = os.path.join(plots_output_dir_path, "coverage_plot.jpeg")
    fig.write_image(image_path)

    return True


def plot_coverage_plots_manager(input_dir, reference_file_path=None):
    if input_dir is None:
        #Take in the reference file as arguments.
        parser = argparse.ArgumentParser()

        #Add argument options
        parser.add_argument("-i", "--input_dir", help="Path to the VarScan consensus table (*.consensus.vcf).", required=True)
        parser.add_argument("-r", "--reference", help="Reference FASTA (used to pad coverage to full amplicon length).", default=None)

        #Fetch the argument data.
        args = vars(parser.parse_args())
        input_dir = args['input_dir']
        reference_file_path = args['reference']

    #Get the consensus file.
    cell_line_id = input_dir.rstrip('/').split('/')[-1]
    # cell_line_id:Spike-Mut10-1_S7_L001.consensus.vcf
    print("[Takeda Pipeline Message] Creating coverage plot for cell line: {} using plotly".format(cell_line_id))

    consensus_file_path = input_dir
    if not os.path.exists(consensus_file_path):
        print("[Takeda Pipeline Message] Sorry, the coverage file is not present. Please make sure a VarScan consensus table (*.consensus.vcf) exists at the path passed to this plot step.")
        return False

    input_dir = '/'.join(input_dir.split('/')[:-1])
    input_dir = os.path.dirname(input_dir)

    #Create a output direcotry for the coverage plot.
    plots_output_dir_path = create_plots_output_directory(input_dir)

    reference_length = None
    if reference_file_path and os.path.exists(reference_file_path):
        reference_length = get_reference_length(reference_file_path)

    consensus_rows = read_consensus_file(consensus_file_path, reference_length)
    plot_graphs(consensus_rows, cell_line_id, plots_output_dir_path)


if __name__ == "__main__":
    plot_coverage_plots_manager(None)
