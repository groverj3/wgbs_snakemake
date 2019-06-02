#!/usr/bin/env python3

# Author: Jeffrey Grover
# Purpose: Calculate the average coverage over an entire genome from mosdepth
# per-base coverage output
# Created: 2/2019

import csv
import gzip
from argparse import ArgumentParser

# Function block


def get_genome_size(fasta_file):
    genome_size = 0
    with open(fasta_file, 'r') as fasta_reader:
        for line in fasta_reader:
            if not line.startswith('>'):
                genome_size += len(line.strip())
    return genome_size


def get_depth(per_base_bedgz):
    total_depth = 0
    n_lines = 0
    with gzip.open(per_base_bedgz, 'rt') as bed_reader:
        for line in bed_reader:
            n_lines += 1
            total_depth += int(line.split('\t')[3])
    average_depth = total_depth / n_lines
    return [total_depth, average_depth]


# Command line Parser

parser = ArgumentParser(
    description='Calculate whole-genome coverage depth in X format from a fasta'
    ' file and the per-base coverage from mosdepth (or similar file).')
parser.add_argument(
    '-f', '--fasta', help='.fasta file for genome', metavar='File')
parser.add_argument(
    '-m', '--mosdepth', help='.bed.gz mosdepth per-base', metavar='File')

fasta_file = parser.parse_args().fasta
per_base_bedgz = parser.parse_args().mosdepth

# Process the files

genome_size = get_genome_size(fasta_file)
depth = get_depth(per_base_bedgz)

coverage = depth[0] / genome_size

# Output to stdout

print('Genome Size:', genome_size, '\t')
print('Total Depth:', depth[0], '\t')
print('Average Depth:', depth[1], '\t')
print('X Coverage:', coverage, '\t')
