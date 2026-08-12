import pandas as pd
import glob
import os
import argparse
import re


def filter_clones_by_position(start_position, end_position=None):
    end_position = end_position or start_position  # if no end_position is provided, consider a single position
    all_txt_files = glob.glob('*.txt')  
    txt_files = [f for f in all_txt_files if not re.search('filtered', f)]  # filter out the files that have 'filtered' in their name

    for txt_file in txt_files:
        with open(txt_file, 'r') as f:
            lines = f.readlines()

        clone_blocks = ''.join(lines).split("CloneID:")
        clone_blocks = clone_blocks[1:]  # The first block is empty string

        if start_position == end_position:
            output_filename = os.path.splitext(txt_file)[0] + f'_filtered_pos_{start_position}.txt'
        else:
            output_filename = os.path.splitext(txt_file)[0] + f'_filtered_pos_{start_position}_{end_position}.txt'

        with open(output_filename, 'w') as out_f:
            for clone_block in clone_blocks:
                block_lines = clone_block.split('\n')

                # Extracting CloneID and SampleIDs
                clone_id = block_lines[0].strip()
                sample_ids = block_lines[1].split(': ')[1].strip()

                # Loading data into a DataFrame
                data = [line.split() for line in block_lines[2:]]
                if len(data) > 1:  # Checking if there are rows apart from the header
                    df = pd.DataFrame(data[1:], columns=data[0])

                    # Converting position to int for comparison
                    df['Position'] = pd.to_numeric(df['Position'], errors='coerce')
                    df = df.dropna(subset=['Position'])
                    df['Position'] = df['Position'].astype(int)

                    # Filtering the DataFrame based on the provided position(s)
                    filtered_df = df[(df['Position'] >= start_position) & (df['Position'] <= end_position)]

                    out_f.write(f"CloneID: {clone_id}\n")
                    out_f.write(f"Replicate Sample IDs: {sample_ids}\n")
                    if filtered_df.empty:
                        out_f.write("No data at provided position range\n")
                    else:
                        out_f.write(f"{filtered_df.to_string(index=False)}\n")
                else:  # If no rows, write the headers and a message
                    out_f.write(f"CloneID: {clone_id}\n")
                    out_f.write(f"Replicate Sample IDs: {sample_ids}\n")
                    out_f.write("No data in this block\n")

                out_f.write("\n")  # To separate blocks



# Setting up command line argument parsing
parser = argparse.ArgumentParser(description='Filter clones by position.')
parser.add_argument('start_position', type=int, help='The start of the position range.')
parser.add_argument('end_position', type=int, nargs='?', help='The end of the position range. Optional.')
args = parser.parse_args()

# Call the function with the position or range of positions provided as command line arguments
filter_clones_by_position(args.start_position, args.end_position)
