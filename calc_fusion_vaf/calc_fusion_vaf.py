#!/usr/bin/env python

"""
calc_fusion_vaf.py
==================
Calculates the variant allele fraction (VAF) of a gene fusion
based on Factera output.

Known Issues
------------
- Currently needs a reference genome to be specified. The feature
    of automatically retrieving the wild-type gene sequences from
    the Ensembl REST API hasn't been implemented yet. This is
    partially due to the problem that gene symbols are passed, not
    stable gene IDs. Therefore, Ensembl might return ambiguous
    genes. Another option is that the breakpoint position can be
    considered when unambiguously determining the gene from Ensembl.
- BWA indexing is computationally expensive, so avoiding re-doing it
    for nothing is an important feature. Now, the script raises an
    error if the output reference file already exists, unless the
    --force_overwrite option is specified. This behaviour can be
    improved.
"""

__author__ = 'Bruno Grande'
__contact__ = 'bgrande@sfu.ca'

import argparse
import csv
import sys
import logging
from pprint import pformat
import os
import shutil
import subprocess
import shlex


def main():
    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_name', '-s', default='calc_fusion_vaf',
                        help='Used as prefix for all generated files in the output directory.')
    parser.add_argument('--gene_1', '-g1', type=str.lower,
                        help='First gene involved in the fusion.')
    parser.add_argument('--gene_2', '-g2', type=str.lower,
                        help='Second gene involved in the fusion.')
    parser.add_argument('--fastq', '-i', nargs=2,
                        help='Location of both FASTQ files (forward and reverse reads).')
    parser.add_argument('--factera_fusions', '-f', type=argparse.FileType('r'),
                        help='Detailed Factera fusions output file (*.fusions.txt).')
    parser.add_argument('--output_dir', '-o',
                        help='Output directory for generated files.')
    parser.add_argument('--reference_fasta', '-r',
                        help=('Location of the reference genome FASTA file. If specified, '
                              'the reference will be included in the alignment in order to '
                              'reduce off-target alignments.'))
    parser.add_argument('--bwa_dir', '-b', default='',
                        help=('Directory in which the BWA binary is located. If not '
                              'specified, this script will assume it\'s in the '
                              'PATH environment variable.'))
    parser.add_argument('--samtools_dir', '-t', default='',
                        help=('Directory in which the samtools binary is located. If not '
                              'specified, this script will assume it\'s in the '
                              'PATH environment variable.'))
    parser.add_argument('--force_overwrite', '-fo', action='store_true',
                        help=('Specify this option in order to overwrite any existing '
                              'output files.'))
    parser.add_argument('--log', '-l', type=str.upper, default='INFO',
                        help='Enable debugging mode.')
    args = parser.parse_args()
    defaultlevel = getattr(logging, 'INFO')
    loglevel = getattr(logging, args.log, defaultlevel)
    logging.basicConfig(
        format='%(levelname)s: %(message)s', stream=sys.stderr, level=loglevel)
    logging.info('Initializing script...')
    logging.debug('Value of args:\n%s', pformat(vars(args)))

    # Parse Factera output and locate relevant fusions
    logging.info('Extracting fusions of specified gene pair...')
    factera_fusions_reader = csv.DictReader(
        args.factera_fusions, delimiter='\t')
    g1_g2_fusions = dict()
    for index, row_dict in enumerate(factera_fusions_reader):
        if ((row_dict['Region1'].lower() == args.gene_1 and
                row_dict['Region2'].lower() == args.gene_2) or
            (row_dict['Region1'].lower() == args.gene_2 and
                row_dict['Region2'].lower() == args.gene_1)):
            row_dict['index'] = index
            g1_g2_fusions[index] = row_dict
    logging.debug('Value of g1_g2_fusions:\n%s', pformat(g1_g2_fusions))

    # Extract sequences for fusions from Factera output
    logging.info('Extracting fusion sequences...')
    # Removes trailing .factera.fusions.txt
    factera_outout_prefix = args.factera_fusions.name[:-20]
    factera_fusionseqs = open(
        factera_outout_prefix + '.factera.fusionseqs.fa', 'r')
    fusion_indices = g1_g2_fusions.keys()
    for seq_dict in fasta_generator(factera_fusionseqs):
        if seq_dict['index'] in fusion_indices:
            g1_g2_fusions[seq_dict['index']]['fusion_name'] = seq_dict['name']
            g1_g2_fusions[seq_dict['index']]['fusion_seq'] = seq_dict['seq']
    logging.debug('Value of g1_g2_fusions:\n%s', pformat(g1_g2_fusions))

    # Generate reference FASTA file for alignment
    logging.info('Generating new reference FASTA file...')
    ## Create output directory if doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    ## Check if reference already exists from previous run
    ## to prevent overwrite
    new_reference_name = args.output_dir + '/' + args.sample_name + '.new_reference.fa'
    if os.path.exists(new_reference_name):
        logging.warning('Output reference FASTA file already exists and presumably other '
                        'output files as well. If --force_overwrite option enabled, these '
                        'will be overwritten. If not, the script will terminate now.')
        if not args.force_overwrite:
            raise IOError('Output file(s) already exist.')
    ## If the original reference genome isn't specified, create new reference
    ## based on the fusion sequences and the wild-type gene sequences
    ## ** Not implemented yet **
    pass
    ## If the original reference genome is specified, append the fusion
    ## sequences to it in a new reference FASTA file
    logging.info('Copying specified reference FASTA file to output directory...')
    shutil.copyfile(args.reference_fasta, new_reference_name)
    logging.info('Appending fusion sequences to new reference FASTA file...')
    new_reference_out = open(new_reference_name, 'a')
    for index, fusion in g1_g2_fusions.items():
        fusion_fasta = '>{Region1}_{Region2}_{Break1}_{Break2}\n{fusion_seq}\n'.format(**fusion)
        new_reference_out.write(fusion_fasta)

    # Run BWA MEM alignment against the new reference
    ## First, create a BWA index for the new reference
    logging.info('Creating BWA index for new reference genome...')
    index_command = '{bwa_dir}/bwa index {new_reference_name}'.format(
        new_reference_name=new_reference_name, **vars(args)
    )
    index_process = subprocess.Popen(shlex.split(index_command))
    index_process.wait()
    ## Second, align FASTQ files to new reference genome
    logging.info('Aligning FASTQ files to new reference genome...')
    fastq_files = ' '.join(args.fastq)
    align_command = '{bwa_dir}/bwa mem {new_reference_name} {fastq_files}'.format(
        new_reference_name=new_reference_name, fastq_files=fastq_files, **vars(args)
    )
    output_bam_name = args.output_dir + '/' + args.sample_name + '.bam'
    sam_to_bam_command = 'samtools view -bT -o {output_bam_name} {new_reference_name} -'.format(
        new_reference_name=new_reference_name, output_bam_name=output_bam_name
    )
    align_process = subprocess.Popen(align_command, stdout=subprocess.PIPE)
    sam_to_bam_process = subprocess.Popen(sam_to_bam_command, stdin=align_process.stdout)
    sam_to_bam_process.wait()


def fasta_generator(file_object):
    """
    Generator for FASTA files, outputting a dictionary
    for each sequence
    """
    index = 0
    seq_dict = None
    # Remove trailing newline
    for line in (raw_line.rstrip('\n') for raw_line in file_object):
        if line.startswith('>'):
            if seq_dict is not None:
                yield seq_dict
                index += 1
            seq_dict = {
                'name': line[1:],
                'index': index,
                'seq': ''
            }
            continue
        seq_dict['seq'] += line
    yield seq_dict


if __name__ == '__main__':
    main()
