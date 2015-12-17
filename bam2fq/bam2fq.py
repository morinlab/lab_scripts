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
-Prints total number of reads processed

Known Issues
-Ignored read shuffling

"""
import os
import argparse
import string
import subprocess
import sys
import glob
import re

def main():
    """Parse BAM file and convert to FASTQ"""
    # Argument parsing
    args = parse_args()
    chnk_size = args.num_reads
    bam_infile = args.bam
    outdir = args.outdir

    # Initialize globals
    paired = False
    p_outstring = '@{0}/1\n{1}\n+\n{2}\n@{0}/2\n{3}\n+\n{4}\n'
    up_outstring = '@{0}\n{1}\n+\n{2}\n'
    bits = [16, 64, 128, 2048]

    chnk_num = 1
    p_read_count = 0
    up_read_count = 0
    tot_read_count = 0
    
    # Initalize data structures
    paired_dict = {}

    # Create output directory if !exist
    if not os.path.exists(os.path.abspath(outdir)):
        os.mkdir(os.path.abspath(outdir))

    outdir = os.path.abspath(outdir)

    gzip_cmd = "gzip > "
    # Open directories for output
    p_sp = spawn_gzip(chnk_num, outdir, gzip_cmd, "paired_")
    chnk_num += 1
    up_sp = spawn_gzip(chnk_num, outdir, gzip_cmd, "unpaired_")

    # Look over each line in bam
    for line in bam_infile:
        bam_line = line.split("\t")
        qname = bam_line[0]
        sam_flag_list = decompose_flag(bam_line[1], bits)
        seq = bam_line[9]
        quality = bam_line[10]

        # Supplementary alignment check
        if sam_flag_list[3]:
            tot_read_count += 1
            continue

        # Reverse complement sequence and reverse quality string
        # if revcomp flag is set
        if sam_flag_list[0]:
            seq = revcomp(seq)
            quality = quality[::-1]

        # Add qname as key if it doesn't exist.
        # If it already exists then current read is part of
        # a pair.
        
        if qname in paired_dict:
            paired = True
        else:
            paired_dict[qname] = ([], [])

        # Add sequence and quality to list corresponding to the
        # read
        if sam_flag_list[1]:
            sort_reads(paired_dict[qname], 0, seq, quality)
        elif sam_flag_list[2]:
            sort_reads(paired_dict[qname], 1, seq, quality)

        if paired:
            t = check_read_count(p_read_count, chnk_num, outdir,
                                 chnk_size, p_sp, gzip_cmd, "paired_")
            p_read_count, chnk_num, p_sp = t
            val = [qname, paired_dict[qname][0][0], paired_dict[qname][0][1],
                   paired_dict[qname][1][0], paired_dict[qname][1][1]]
            p_sp.stdin.write(p_outstring.format(*val))
            p_read_count += 2
            tot_read_count += 2
            del paired_dict[qname]
            paired = False
    else:
        bam_infile.close()

    # Output unpaired reads
    for qname in paired_dict.keys():
        t = check_read_count(up_read_count, chnk_num, outdir,
                             chnk_size, up_sp, gzip_cmd, "unpaired_")
        up_read_count , chnk_num, up_sp = t
        val = [qname]

        if len(paired_dict[qname][0]):
            val = val + paired_dict[qname][0]
        else:
            val = val + paired_dict[qname][1]

        up_sp.stdin.write(up_outstring.format(*val))
        up_read_count += 1
        tot_read_count += 1
        del paired_dict[qname]

    # Write interval file
    write_interval_file(outdir)
    return


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


def spawn_gzip(chnk_num, chnk_dir, gzip_cmd, read_type):
    """Returns gzip process filehandle."""
    s = os.path.join(chnk_dir, str(read_type) + "chunk" + str(chnk_num) + '.fastq.gz')
    p = subprocess.Popen(gzip_cmd + s, bufsize=-1, shell=True,
                         stdin=subprocess.PIPE)
    return p


def write_interval_file(directory):
    """Creates an interval file in the given directory"""
    interval_filename = os.path.join(directory, 'interval.txt')
    files = glob.glob(os.path.join(directory, "*.fastq.gz"))
    outstring = [str(x)[:-len(".fastq.gz")] + '\n' for x in files]
    interval_file = open(interval_filename, 'w')
    interval_file.write(string.join(outstring, ''))
    interval_file.close()
    return


def check_read_count(rc, chnk_num, chnk_dir, chnk_size, sp, cmd, read_type):
    """Checks if read count limit reached, spawns new gzip process."""
    if rc >= chnk_size:
        chnk_num += 1
        sp = spawn_gzip(chnk_num, chnk_dir, cmd, read_type)
        rc = 0
    return [rc, chnk_num, sp]


def parse_args():
    """Parse command-line arguments"""
    # Setup
    parser = argparse.ArgumentParser()
    # Optional arguments
    parser.add_argument('--num_reads', '-n', type=int, default=75000000,
                        help='Specify number of reads per FASTQ file.')
    # Positional arguments
    parser.add_argument('bam', type=argparse.FileType('r'),
                        help='Specify SAM file.')
    parser.add_argument("outdir", default='.')
    # Parsing
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
