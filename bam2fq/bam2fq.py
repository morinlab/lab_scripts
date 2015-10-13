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
-Ignored read shuffling

"""

import argparse
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
    
    # Initialize globals
    paired = False
    p_outstring = '@{0}/1\n{1}\n+\n{2}\n@{0}/2\n{3}\n+\n{4}\n'
    up_outstring = '@{0}\n{1}\n+\n{2}\n'
    bits = [1, 16, 64, 128, 256, 512, 1024, 2048]
    p_chunk_number = 1
    p_read_count = 0
    up_chunk_number = 1
    up_read_count = 0
    
    # Initalize data structures
    paired_dict = {}

    # Creates output directories if they don't exist
    if not os.path.exists(unpaired_outdir):
        os.mkdir(unpaired_outdir)

    if not os.path.exists(paired_outdir):
        os.mkdir(paired_outdir)

    # Open directories for output
    paired_chunk_file = open_chunk_file(p_chunk_number, paired_outdir)
    unpaired_chunk_file = open_chunk_file(up_chunk_number, unpaired_outdir)
    
    # Look over each line in bam
    for line in bam_infile:
        bam_line = line.split("\t")
        qname = bam_line[0]
        sam_flag_list = decompose_flag(bam_line[1], bits)
        seq = bam_line[9]
        quality = bam_line[10]

        # Reverse complement sequence and reverse quality string
        # if revcomp flag is set
        if sam_flag_list[1]:
            seq = revcomp(seq)
            quality = quality[::-1]

        # Unpaired read if supplementary flag is set    
        if sam_flag_list[4] or sam_flag_list[5] or sam_flag_list[6] or sam_flag_list[7]:
            if up_read_count >= chunk_size:
                unpaired_chunk_file.close()
                up_chunk_number += 1
                up_read_count = 0
                unpaired_chunk_file = open_chunk_file(up_chunk_number, unpaired_outdir)
            read = [qname, seq, quality]
            write_to_chunk(unpaired_chunk_file, read, up_outstring)
            up_read_count += 1
            continue

        # Add qname as key if it doesn't exist.
        # If it already exists then current read is part of
        # a pair.
        if qname not in paired_dict.keys():
            paired_dict[qname] = ([], [])
        else:
            paired = True

        # Add sequence and quality to list corresponding to the
        # read
        if sam_flag_list[2]:
            sort_reads(paired_dict[qname], 0, seq, quality)
        elif sam_flag_list[3]:
            sort_reads(paired_dict[qname], 1, seq, quality)

        if paired:       
            if p_read_count >= chunk_size:
                paired_chunk_file.close()
                p_chunk_number += 1
                p_read_count = 0
                paired_chunk_file = open_chunk_file(p_chunk_number, paired_outdir)
            t = [qname, paired_dict[qname][0][0], paired_dict[qname][0][1],
                 paired_dict[qname][1][0], paired_dict[qname][1][1]]
            write_to_chunk(paired_chunk_file, t, p_outstring)
            p_read_count += 2
            del paired_dict[qname]
            paired = False
    else:
        bam_infile.close()
        paired_chunk_file.close()
        create_interval_file(p_chunk_number, paired_outdir, 'paired')
    
    # Output unpaired reads
    for qname in paired_dict.keys():
        if up_read_count >= chunk_size:
            unpaired_chunk_file.close()
            up_chunk_number += 1
            up_read_count = 0
            unpaired_chunk_file = open_chunk_file(up_chunk_number, unpaired_outdir)
        read = [qname]
        if len(paired_dict[qname][0]):
            read = read + paired_dict[qname][0]
        else:
            read = read + paired_dict[qname][1]
        write_to_chunk(unpaired_chunk_file, read, up_outstring)
        up_read_count += 1
        del paired_dict[qname]
    else:
        unpaired_chunk_file.close()
        create_interval_file(up_chunk_number, unpaired_outdir, 'unpaired')


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
    """DEPRECATED"""
    """Iterates to next chunk if max. reads has been reached"""
    flushed = False
    if read_count > chunk_size:
        chunk_file.flush()
        chunk_file.close()
        chunk_number += 1
        read_count = 0
        chunk_file = open_chunk_file(chunk_number, output_dir)
        flushed = True
    return chunk_file, read_count, chunk_number, flushed


def create_interval_file(chunk_number, directory, prefix):
    """Creates an interval file in the given directory"""
    interval_filename = os.path.join(directory, prefix + '_interval.txt')
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
    import timeit
    print(timeit.timeit("main()", setup="from __main__ import main", number=1))
