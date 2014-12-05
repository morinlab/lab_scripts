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
- This was meant for capture-based sequencing data, where one of the
    two genes in a fusion is targeted. This has implications for the
    calculation of the VAF, where the coverage for the wild-type
    alleles needs to be doubled since coverage for one of the genes
    is missing.

Known Issues
------------
- Currently needs a reference genome to be specified. The feature
    of automatically retrieving the wild-type gene sequences from
    the Ensembl REST API hasn't been implemented yet. This is
    partially due to the problem that gene symbols are passed, not
    stable gene IDs. Therefore, Ensembl might return ambiguous
    genes. Another option is that the breakpoint position can be
    considered when unambiguously determining the gene from Ensembl.
    Also, the breakpoint positions are based on the reference.
- Currently, reads that span the edges of the window are included.
    I believe it would be better if the reads need to be entirely
    within the window.
- I do not remove duplicates.
"""

__author__ = 'Bruno Grande'
__contact__ = 'bgrande@sfu.ca'

import argparse
import csv
import sys
import logging
import os
import shutil
import subprocess
import pysam
from collections import defaultdict

OUTPUT_FIELDNAMES = [
    'fusion', 'gene_1', 'gene_2', 'breakpoint_1', 'breakpoint_2', 'fusion_spanning_reads',
    'fusion_spanning_read_pairs', 'fusion_total_support', 'wildtype_spanning_reads',
    'wildtype_spanning_read_pairs', 'wildtype_total_support', 'variant_allele_fraction']


def main():
    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('reference_genome', help='Location of reference genome in FASTA format')
    parser.add_argument('fastq_files', nargs=2, help='Input FASTQ files')
    parser.add_argument('factera_fusions_txt', type=argparse.FileType('r'),
                        help='Detailed Factera fusions output file (*.fusions.txt).')
    parser.add_argument('--window', '-w', type=int, default=2000,
                        help='Size of window around breakpoint positions')
    parser.add_argument('--min_overlap', '-m', type=int, default=10,
                        help='Minimum number of overlapping bases for a spannning read.')
    parser.add_argument('--output_dir', '-o', default='.', help='Output directory')
    parser.add_argument('--log', '-l', help='Specify log file. Otherwise, stderr.')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads used')
    args = parser.parse_args()
    # Setup logging
    if args.log:
        logging.basicConfig(
            format='%(levelname)s: %(message)s', level=logging.INFO, filename=args.log)
    else:
        logging.basicConfig(
            format='%(levelname)s: %(message)s', level=logging.INFO, stream=sys.stderr)
    logging.info('Initializing script...')

    # Check if final output file already exists
    results_file_name = args.output_dir + '/results.txt'
    if os.path.exists(results_file_name):
        logging.error('Final output file already exists. Exiting...')
        sys.exit()

    # Parse Factera output and locate relevant fusions
    logging.info('Extracting fusions of specified gene pair...')
    factera_fusions_reader = csv.DictReader(
        args.factera_fusions_txt, delimiter='\t')
    fusions = dict()
    for index, row_dict in enumerate(factera_fusions_reader):
        row_dict['index'] = index
        fusions[index] = row_dict

    # Extract sequences for fusions from Factera output
    logging.info('Extracting fusion sequences...')
    # Removes trailing .factera.fusions.txt
    factera_outout_prefix = args.factera_fusions_txt.name[:-20]
    factera_fusionseqs = open(
        factera_outout_prefix + '.factera.fusionseqs.fa', 'r')
    for seq_dict in fasta_gen(factera_fusionseqs):
        fusions[seq_dict['index']]['fusion_name'] = seq_dict['name']
        fusions[seq_dict['index']]['fusion_seq'] = seq_dict['seq']
        fusions[seq_dict['index']]['fusion_ref_name'] = (
            '{Region1}_{Region2}_{Break1}_{Break2}'.format(**fusions[seq_dict['index']])
        )

    # Generate reference FASTA file for alignment
    ## Create output directory if doesn't exist
    if not os.path.exists(args.output_dir):
        logging.info('Creating output directory...')
        os.makedirs(args.output_dir)
    else:
        logging.warning('Specified output directory already exists...')
    ## Check if reference already exists from previous run
    ## to prevent overwrite
    new_reference_name = args.output_dir + '/reference_genome.with_fusions.fa'
    if os.path.exists(new_reference_name):
        logging.warning('Reference genome with fusions already exists. Skipping...')
    else:
        logging.info('Generating new reference genome with fusion sequences...')
        ## Append the fusion sequences to it in a new reference FASTA file
        logging.info('Copying specified reference FASTA file to output directory...')
        shutil.copyfile(args.reference_genome, new_reference_name)
        logging.info('Appending fusion sequences to new reference FASTA file...')
        with open(new_reference_name, 'a') as new_reference_out:
            for index, fusion in fusions.items():
                fusion_fasta = '>{fusion_ref_name}\n{fusion_seq}\n'.format(**fusion)
                new_reference_out.write(fusion_fasta)

    # Run BWA MEM alignment against the new reference
    ## Create a BWA index for the new reference
    if (os.path.exists(new_reference_name + '.amb') and
            os.path.exists(new_reference_name + '.ann') and
            os.path.exists(new_reference_name + '.bwt') and
            os.path.exists(new_reference_name + '.pac') and
            os.path.exists(new_reference_name + '.sa')):
        logging.warning('Index for reference genome with fusions already exists. Skipping...')
    else:
        logging.info('Creating BWA index for new reference genome...')
        index_cmd = ['bwa', 'index', new_reference_name]
        run_cmd(index_cmd)
    ## Align FASTQ files to new reference genome
    output_sai_1 = args.output_dir + '/bwa_aln_1.sai'
    output_sai_2 = args.output_dir + '/bwa_aln_2.sai'
    output_sam = args.output_dir + '/bwa_sampe.sam'
    output_bam = args.output_dir + '/bwa_sampe.bam'
    output_bam_sorted = args.output_dir + '/bwa_sampe.sorted.bam'
    ### Running BWA aln
    if os.path.exists(output_sai_1) and os.path.exists(output_sai_2):
        logging.warning('BWA aln output already exists. Skipping...')
    else:
        logging.info('Running BWA aln for both FASTQ files...')
        align_cmd_prefix = ['bwa', 'aln', '-t', args.threads, '-f']
        align_cmd_1 = align_cmd_prefix + [output_sai_1, new_reference_name, args.fastq_files[0]]
        align_cmd_2 = align_cmd_prefix + [output_sai_2, new_reference_name, args.fastq_files[1]]
        run_cmd(align_cmd_1)
        run_cmd(align_cmd_2)
    ### Running BWA sampe
    if os.path.exists(output_sam):
        logging.warning('BWA sampe output already exists. Skipping...')
    else:
        logging.info('Running BWA sampe...')
        sampe_cmd = [
            'bwa', 'sampe', '-f', output_sam, new_reference_name, output_sai_1, output_sai_2,
            args.fastq_files[0], args.fastq_files[1]]
        run_cmd(sampe_cmd)
    ### Converting from SAM to BAM format
    if os.path.exists(output_bam):
        logging.warning('Converted BAM file already exists. Skipping...')
    else:
        logging.info('Converting SAM file to BAM format...')
        pysam.view('-S', '-b', '-o' + output_bam, output_sam)
    ### Sorting converted BAM file
    if os.path.exists(output_bam_sorted):
        logging.warning('Sorted BAM file already exists. Skipping...')
    else:
        logging.info('Sorting converted BAM file...')
        pysam.sort('-f', output_bam, output_bam_sorted)
    ### Indexing sorted BAM file
    if os.path.exists(output_bam_sorted + '.bai'):
        logging.warning('Index for sorted BAM file already exists. Skipping...')
    else:
        logging.info('Indexing sorted output BAM file...')
        pysam.index(output_bam_sorted)

    # Quantify support for fusion alleles
    logging.info('Loading generated BAM file for analysis...')
    bam_file = pysam.AlignmentFile(output_bam_sorted, 'rb')
    ## Extract chromosome lengths in order to calculate fusion middle point
    chr_lengths = dict(zip(bam_file.references, bam_file.lengths))
    ## Calculate coverage values for each fusion and wild-type alleles
    logging.info('Calculating coverage for fusion and wild-type sequences...')
    results = []
    for index, fusion in fusions.items():
        # Calculating coverage for the fusion
        position = chr_lengths[fusion['fusion_ref_name']] / 2
        spanning_reads, spanning_read_pairs = calc_cov_at_pos(
            bam_file, fusion['fusion_ref_name'], position, args.window, args.min_overlap)
        fusion_spanning_reads = len(spanning_reads)
        fusion_spanning_read_pairs = len(spanning_read_pairs)
        # Calculating coverage for the wild-type allele
        bp1_chr, bp1_pos = fusion['Break1'].split(':')
        bp2_chr, bp2_pos = fusion['Break2'].split(':')
        bp1_cov = calc_cov_at_pos(bam_file, bp1_chr, bp1_pos, args.window, args.min_overlap)
        bp2_cov = calc_cov_at_pos(bam_file, bp2_chr, bp2_pos, args.window, args.min_overlap)
        wildtype_spanning_reads = bp1_cov[0] + bp2_cov[0]
        wildtype_spanning_read_pairs = bp1_cov[1] + bp2_cov[1]
        # Calculate support and vaf
        fusion_total_support = fusion_spanning_reads + fusion_spanning_read_pairs
        wildtype_total_support = wildtype_spanning_reads + wildtype_spanning_read_pairs
        variant_allele_fraction = round(
            float(fusion_total_support) / float(wildtype_total_support), 1)
        # Compile results
        results.append({
            'fusion': fusion['fusion_ref_name'],
            'gene_1': fusion['Region1'],
            'gene_2': fusion['Region2'],
            'breakpoint_1': fusion['Break1'],
            'breakpoint_2': fusion['Break2'],
            'fusion_spanning_reads': fusion_spanning_reads,
            'fusion_spanning_read_pairs': fusion_spanning_read_pairs,
            'fusion_total_support': fusion_total_support,
            'wildtype_spanning_reads': wildtype_spanning_reads,
            'wildtype_spanning_read_pairs': wildtype_spanning_read_pairs,
            'wildtype_total_support': wildtype_total_support,
            'variant_allele_fraction': variant_allele_fraction
        })

    # Write out results to output file
    logging.info('Outputting results...')
    with open(results_file_name, 'w') as results_file:
        csv_writer = csv.DictWriter(results_file, OUTPUT_FIELDNAMES)
        csv_writer.writerows(results)

    # Clean up
    os.remove(output_sai_1)
    os.remove(output_sai_2)
    os.remove(output_sam)
    os.remove(output_bam)


def fasta_gen(file_object):
    """
    Generator for FASTA files, outputting a dictionary
    for each sequence.
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
    Standardizes the way commands are run.
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


def calc_cov_at_pos(aln_file, chromosome, pos, window=2000, min_overlap=10):
    """
    Calculates the number of reads and read pairs that span
    a specified position within a given window.
    Returns tuple of lists, (spanning_reads, spanning_read_pairs):
    - spanning_reads contains spanning reads
    - spanning_read_pairs contains tuples of read pairs
    """
    spanning_reads = []
    spanning_read_pairs = []
    non_spanning_reads = defaultdict(list)
    # Extract reads from specified region from alignment file
    start = pos - window / 2
    end = pos + window / 2
    for read in aln_file.fetch(chromosome, start, end):
        # Check if read spans breakpoint
        if (read.reference_start < pos - min_overlap and read.reference_end > pos + min_overlap):
            spanning_reads.append(read)
        else:
            # If not, cache non-spanning reads in dict according to read name
            non_spanning_reads[read.query_name].append(read)
    # Go through cached reads while only considered pairs
    for r1, r2 in ((x[0], x[1]) for x in non_spanning_reads.values() if len(x) == 2):
        # Check if reads are on opposite strands
        if (r1.is_reverse and not r2.is_reverse) or (not r1.is_reverse and r2.is_reverse):
            # Figure out which read is on positive strand and vice versa
            if r1.is_reverse:
                r_plus, r_minus = r2, r1
            else:
                r_plus, r_minus = r1, r2
            # Check if pointing each other, i.e. positive insert size
            if r_plus.template_length > 0:
                # Check if breakpoint is between both reads, i.e. not spanning
                if ((r_plus.reference_start + 1 < pos - min_overlap) and
                        (r_minus.reference_start + 1 + r_minus.reference_length >
                         pos + min_overlap)):
                    spanning_read_pairs.append((r1, r2))
                    continue
    return (spanning_reads, spanning_read_pairs)


if __name__ == '__main__':
    main()
