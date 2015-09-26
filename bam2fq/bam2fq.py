#!/usr/bin/env python

"""
bam2fq.py
==========
Description
Converts each read in a SAM file into FASTQ format while pairing
together mate end pair reads. It will output multiple FASTQ files
each containing a maximum number of reads (default: 75,000,000).

Dependencies
-None

Inputs
-SAM-formatted file of read mappings
-Output directory for FASTQ chunks

Outputs
-FASTQ files for paired reads
-FASTQ files for unpaired reads

Known Issues
-BUG: if the user sets --num_reads to an odd number (e.g., 5), 
but there are an even number of read pairs (e.g., 3 pairs or 6 reads), 
the algorithm will keep the read pairs together rather than splitting 
them into separate files (i.e., will create one file with 6 reads 
rather than one file with 5 reads and another file with 1 read)

"""

import argparse
import random
import os
import gzip
import string


def main():
    """Parse BAM file and convert to FASTQ"""
    # Argument parsing
    args = parse_args()
    chunk_size = args.num_reads
    bam_infile = args.bam[0]
    unpaired_outdir = args.single_outdir[0]
    paired_outdir = args.paired_outdir[0]

    # Initalize data structures
    read_pair_dict = {}
    paired_read_names = []

    # Creates output directories if they don't exist
    if not os.path.exists(unpaired_outdir):
        os.mkdir(unpaired_outdir)

    if not os.path.exists(paired_outdir):
        os.mkdir(paired_outdir)

    # Look over each line in bam
    for line in bam_infile:
        bam_line = line.split("\t")
        qname = bam_line[0]
        sam_flag_list = decompose_flag(bam_line[1], [1, 16, 64, 128])
        seq = bam_line[9]
        quality = bam_line[10]

        # Reverse complement sequence and reverse quality string
        # if revcomp flag is set
        if sam_flag_list[1]:
            seq = revcomp(seq)
            quality = quality[::-1]

        # Add qname as key if it doesn't exist.
        # If it already exists, means the current read is part of
        # a pair.
        if qname not in read_pair_dict.keys():
            read_pair_dict[qname] = ([], [])
        else:
            paired_read_names.append(qname)

        # Add sequence and quality to array corresponding to the
        # read
        if sam_flag_list[2]:
            sort_reads(read_pair_dict[qname], 0, seq, quality)
        elif sam_flag_list[3]:
            sort_reads(read_pair_dict[qname], 1, seq, quality)
    else:
        bam_infile.close()

    if len(paired_read_names):
        chunk_number = 1
        read_count = 0
        paired_chunk_file = open_chunk_file(chunk_number, paired_outdir)
        outstring = '@{0}/1\n{1}\n+\n{2}\n@{0}/2\n{3}\n+\n{4}\n'
        random.shuffle(paired_read_names)
        for k in paired_read_names:
            t = [k, read_pair_dict[k][0][0], read_pair_dict[k][0][1],
                 read_pair_dict[k][1][0], read_pair_dict[k][1][1]]
            write_to_chunk(paired_chunk_file, t, outstring)
            read_count += 2
            paired_chunk_file, read_count, chunk_number =\
            check_read_count(paired_chunk_file, read_count, chunk_number,
                             chunk_size, paired_outdir)
        # Take current chunk_number and format a string 
        create_interval_file(chunk_number, paired_outdir, 'paired')

    # Remaining keys in dictionary should be unpaired reads
    unpaired_read_names = set(read_pair_dict.keys()) - set(paired_read_names)
    if len(unpaired_read_names):
        chunk_number = 1
        read_count = 0
        unpaired_chunk_file = open_chunk_file(chunk_number, unpaired_outdir)
        outstring = '@{0}\n{1}\n+\n{2}'
        for k in unpaired_read_names:
            read = [k]

            # Determine which read (fwd or rev) is in the dictionary
            if len(read_pair_dict[k][0]):
                read = read + read_pair_dict[k][0]
            else:
                read = read + read_pair_dict[k][1]

            write_to_chunk(unpaired_chunk_file, read, outstring)
            read_count += 1
            unpaired_chunk_file, read_count, chunk_number =\
            check_read_count(unpaired_chunk_file, read_count,
                             chunk_number, chunk_size, unpaired_outdir)
        create_interval_file(chunk_number, unpaired_outdir, 'unpaired')

def decompose_flag(flag, bits=[2**x for x in range(12)]):
    """Identifies the bits that are set in the SAM flag"""
    flag_list = [bool(int(flag) & x) for x in bits]
    return flag_list


def revcomp(seq):
    """Reverse complements a sequence"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rev = "".join(complement.get(base, base) for base in reversed(seq))
    return rev


def sort_reads(tuple_to_append, list_index, seq, quality):
    """Appends items to a list at 'list_index' in a given tuple reference"""
    tuple_to_append[list_index].append(seq)
    tuple_to_append[list_index].append(quality)


def open_chunk_file(chunk_number, filepath="."):
    """Opens new chunk file for writing"""
    chunk_file = os.path.join(filepath, '{0}.fastq.gz'.format(chunk_number))
    return gzip.open(chunk_file, 'wb')


def write_to_chunk(chunk_file, read_data, outstring):
	"""Writes formatted string to current chunk"""
	chunk_file.write(outstring.format(*read_data))


def check_read_count(chunk_file, read_count, chunk_number, chunk_size,
                     output_dir):
    """Iterates to next chunk if max. reads has been reached"""
    if read_count >= chunk_size:
		chunk_file.flush()
		chunk_file.close()
		chunk_number += 1
		read_count = 0
		chunk_file = open_chunk_file(chunk_number, output_dir)
    return chunk_file, read_count, chunk_number


def create_interval_file(chunk_number, directory, prefix):
    """Creates an interval file in the given directory"""
    interval_filename = os.path.join(directory,
                                     '{0}_interval.txt'.format(prefix))
    outstring = string.join([str(i) + '.fastq.gz' for i in range(1, chunk_number + 1)], '\n')
    interval_file = open(interval_filename, 'w')
    interval_file.write(outstring)
    interval_file.close()


def parse_args():
    """Parse command-line arguments"""
    # Setup
    parser = argparse.ArgumentParser()
    # Optional arguments
    parser.add_argument('--num_reads', '-n', type=int, default=75000000,
                        help='Specify number of reads per FASTQ file.')
    # Positional arguments
    parser.add_argument('bam', nargs=1, type=argparse.FileType('r'),
                        help='Specify SAM file.')
    parser.add_argument("paired_outdir", default='.', nargs=1)
    parser.add_argument("single_outdir", default='.', nargs=1)
    # Parsing
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
