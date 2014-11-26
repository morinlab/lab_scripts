#!/usr/bin/env python


import argparse


def main():
    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_name', '-s',
                        help='Used as prefix for all generated files in the output directory.')
    parser.add_argument('--gene_1', '-g1',
                        help='First gene involved in the fusion.')
    parser.add_argument('--gene_2', '-g2',
                        help='Second gene involved in the fusion.')
    parser.add_argument('--fastq', '-i', nargs=2
                        help='Location of both FASTQ files (forward and reverse reads).')
    parser.add_argument('--factera_output', '-f',
                        help='Detailed Factera output file (*.fusions.txt).')
    parser.add_argument('--output_dir', '-o',
                        help='Output directory for generated files.')
    parser.add_argument('--reference', '-r',
                        help=('Location of the reference genome FASTA file. If specified, '
                              'the reference will be included in the alignment in order to '
                              'reduce off-target alignments.'))
    args = parser.parse_args()


if __name__ == '__main__':
    main()
