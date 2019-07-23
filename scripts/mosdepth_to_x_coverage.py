#!/usr/bin/env python3

# Author: Jeffrey Grover
# Purpose: Calculate the average coverage over an entire genome from mosdepth
# per-base coverage output
# Created: 2019-07-10

import gzip
from argparse import ArgumentParser


# Functions

def get_genome_size(fasta_file):
    genome_size = 0
    with open(fasta_file, 'r') as fasta_reader:
        for line in fasta_reader:
            if not line.startswith('>'):
                genome_size += len(line.strip())
    return genome_size


def get_depth(coverage_bedgz):
    bases_sequenced = 0
    with gzip.open(coverage_bedgz, 'rt') as bed_reader:
        for line in bed_reader:
            length = int(line.split('\t')[2]) - int(line.split('\t')[1])
            bases_sequenced += int(line.split('\t')[3]) * length
    return bases_sequenced


def get_x_coverage(bases_sequenced, genome_size):
    return bases_sequenced / genome_size


# Command line Parser

def get_args():
    parser = ArgumentParser(
        description='Calculate X coverage from mosdepth .bed.gz output. Works '
        'over window-based output or per-base output.')
    parser.add_argument(
        '-f', '--fasta',
        help='.fasta file for genome',
        metavar='File')
    parser.add_argument(
        '-m', '--mosdepth',
        help='.bed.gz mosdepth output',
        metavar='File')
    return parser.parse_args()


# Process the files

def main(args):
    genome_size = get_genome_size(args.fasta)
    depth = get_depth(args.mosdepth)
    x_coverage = get_x_coverage(depth, genome_size)

    # Output to stdout

    print('Genome Size:', genome_size)
    print('Total Depth:', depth)
    print('X Coverage:', x_coverage)


if __name__ == '__main__':
    args = get_args()
    main(args)
