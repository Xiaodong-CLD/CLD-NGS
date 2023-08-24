#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filepath', help='Path of coverage file')
parser.add_argument('-o', help='Output filepath', default='coverage_1.txt')
args = parser.parse_args()

with open(args.filepath, 'r') as fi, open(args.o, 'w') as fo:
    fi.readline()
    fo.write('name\tcoverage\tmeandepth\n') 
    for line in fi:
        line = line.rstrip().split('\t')
        coverage = line[5]
        if float(coverage) > 1:
            virus_name = line[0]
            depth = line[6]
            fo.write(f"{virus_name}\t{coverage}\t{depth}\n")

