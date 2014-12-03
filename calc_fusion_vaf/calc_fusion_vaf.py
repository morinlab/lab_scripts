#!/usr/bin/env python

"""
calc_fusion_vaf.py
==================
Calculates the variant allele fraction (VAF) of a gene fusion
based on Factera output.

Assumptions
-----------
- The breakpoint in the fusion sequences is right in the middle,
    which is what Factera does by default.
- BWA and samtools should be located in the PATH environment variable.

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
import pysam


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
    parser.add_argument('--window', '-w', type=int, default=1000,
                        help=('Number of bases on each side of the breakpoint from which'
                              'reads are extracted.'))
    parser.add_argument('--min_overlap', '-mo', type=int, default=10,
                        help=('Minimum number of overlapping bases for a spannning '
                              'read.'))
    parser.add_argument('--output_dir', '-o',
                        help='Output directory for generated files.')
    parser.add_argument('--reference_fasta', '-r',
                        help=('Location of the reference genome FASTA file. If specified, '
                              'the reference will be included in the alignment in order to '
                              'reduce off-target alignments.'))
    parser.add_argument('--force_overwrite', '-fo', action='store_true',
                        help=('Specify this option in order to overwrite any existing '
                              'output files.'))
    parser.add_argument('--log', '-l', help='Specify log file. Otherwise, stderr.')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads used.')
    args = parser.parse_args()
    if args.log:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO,
                            filename=args.log)
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO,
                            stream=sys.stderr)
    logging.info('Initializing script...')

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
    with open(new_reference_name, 'a') as new_reference_out:
        for index, fusion in g1_g2_fusions.items():
            fusion['fusion_ref_name'] = '{Region1}_{Region2}_{Break1}_{Break2}'.format(**fusion)
            fusion_fasta = '>{fusion_ref_name}\n{fusion_seq}\n'.format(**fusion)
            new_reference_out.write(fusion_fasta)

    # Run BWA MEM alignment against the new reference
    ## First, create a BWA index for the new reference
    logging.info('Creating BWA index for new reference genome...')
    index_cmd = ['bwa', 'index', new_reference_name]
    run_cmd(index_cmd)
    ## Second, align FASTQ files to new reference genome
    ### Running BWA aln
    logging.info('Aligning FASTQ files to new reference genome...')
    align_cmd_prefix = ['bwa', 'aln', '-t', args.threads, '-f']
    output_sai_1 = args.output_dir + '/' + args.sample_name + '.1.sai'
    output_sai_2 = args.output_dir + '/' + args.sample_name + '.2.sai'
    align_cmd_1 = align_cmd_prefix + [output_sai_1, new_reference_name, args.fastq[0]]
    align_cmd_2 = align_cmd_prefix + [output_sai_2, new_reference_name, args.fastq[1]]
    run_cmd(align_cmd_1)
    run_cmd(align_cmd_2)
    ### Running BWA sampe
    output_sam = args.output_dir + '/' + args.sample_name + '.sam'
    output_bam = args.output_dir + '/' + args.sample_name + '.bam'
    sampe_cmd = [
        'bwa', 'sampe', '-f', output_sam, new_reference_name, output_sai_1, output_sai_2,
        args.fastq[0], args.fastq[1]
    ]
    run_cmd(sampe_cmd)
    logging.info('Converting SAM file to BAM format...')
    pysam.view('-S', '-b', '-o' + output_bam, output_sam)
    logging.info('Sorting output BAM file...')
    output_bam_sorted = args.output_dir + '/' + args.sample_name + '.sorted.bam'
    pysam.sort('-f', output_bam, output_bam_sorted)
    logging.info('Indexing sorted output BAM file...')
    pysam.index(output_bam_sorted)

    # Quantify support for fusion breakpoint
    ## Iterate through each of the fusion references ("pseudo-chromosomes")
    logging.info('Loading generated BAM file for analysis...')
    output_bam_sorted = ('/Users/bgrande/Desktop/calc_fusion_vaf_test/'
                         'ewings_sarcoma_test.sorted.bam')
    bam_file = pysam.AlignmentFile(output_bam_sorted, 'rb')
    fusion_lengths = dict(zip(bam_file.references, bam_file.lengths))
    print fusion_lengths
    total_spanning_reads = []
    total_spanning_read_pairs = []
    # for index, fusion in g1_g2_fusions.items():
    #     # Extract reads in target window
    #     middle = fusion_lengths[fusion['fusion_ref_name']] / 2
    #     start = middle - args.window
    #     end = middle + args.window
    #     for read in bam_file.fetch(fusion['fusion_ref_name'], start, end):
    #         # Check if read spans breakpoint
    #         if (read.reference_start < middle - args.min_overlap and
    #                 read.reference_end > middle + args.min_overlap):
    #             total_spanning_reads.append(read)
    #     print len(total_spanning_reads)
    bam_file.close()


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


def run_cmd(cmd_args):
    """
    Standardizes the way commands are run
    """
    # Ensure that all args are strings
    cmd_args = [str(x) for x in cmd_args]
    # Log the command being run
    logging.info('Running command:\n {}'.format(' '.join(cmd_args)))
    # Run the command
    cmd_proc = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Capture stdout and stderr, and output to log
    stdout, stderr = cmd_proc.communicate()
    logging.info(stderr)


if __name__ == '__main__':
    main()
