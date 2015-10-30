#!/usr/bin/env python

"""
simple_bam2fq.py
==========
Description
Converts each read in a name sorted SAM file into FASTQ format while
pairing together mate end pair reads. It will output multiple FASTQ
files each containing a maximum number of reads (default: 5,000,000,000).

Dependencies
-None

Inputs
-Name sorted SAM-formatted file of read mappings
-Output directory for FASTQ chunks

Outputs
-FASTQ files for paired reads
-FASTQ files for unpaired reads

Known Issues
-None

"""

import argparse
import os
import subprocess
import string


def main():

    args = parse_args()
    chnk_size = args.num_reads
    bam_infile = args.bam[0]
    outdir = args.outdir[0]

    p_outstring = '@{0}/1\n{1}\n+\n{2}\n@{0}/2\n{3}\n+\n{4}\n'
    up_outstring = '@{0}\n{1}\n+\n{2}\n'

    p_read_count = 0
    up_read_count = 0
    chnk_num = 1

    tot_read_count = 0

    bits = [16, 64, 128, 2048]
    if not os.path.exists(os.path.abspath(outdir)):
        os.mkdir(os.path.abspath(outdir))

    outdir = os.path.abspath(outdir)

    gzip_cmd = "gzip > "
    p_sp = spawn_gzip(chnk_num, outdir, gzip_cmd, "paired_")
    chnk_num += 1
    up_sp = spawn_gzip(chnk_num, outdir, gzip_cmd, "unpaired_")
    # prev_line = [qname, seq, quality, bit 16, bit 64, bit 128]
    prev_line = []
    for line in bam_infile:
        bam_line = line.split('\t')
        sam_flag_list = decompose_flag(bam_line[1], bits)

        # Supplementary alignment check
        if sam_flag_list[3]:
            tot_read_count += 1
            continue

        if sam_flag_list[0]:
            bam_line[9] = revcomp(bam_line[9])
            bam_line[10] = str(bam_line[10])[::-1]

        if not len(prev_line):
            prev_line = [bam_line[0]] + bam_line[9:11] + sam_flag_list[0:3]
            continue

        if str(bam_line[0]) == str(prev_line[0]):
            # check readcount, get process, write to its stdin
            t = check_read_count(p_read_count, chnk_num, outdir,
                                 chnk_size, p_sp, gzip_cmd, "paired_")
            p_read_count, chnk_num, p_sp = t
            val = []
            if sam_flag_list[1]:
                # curr read is first in pair
                val = [bam_line[0], bam_line[9], bam_line[10],
                       prev_line[1], prev_line[2]]
            else:
                # curr read is second in pair
                val = [prev_line[0], prev_line[1], prev_line[2],
                       bam_line[9], bam_line[10]]
            p_sp.stdin.write(p_outstring.format(*val))
            p_read_count += 2
            tot_read_count += 2
            prev_line = []
        else:
            # save prev_line as unpaired read, curr line becomes prev_line
            t = check_read_count(up_read_count, chnk_num, outdir,
                                 chnk_size, up_sp, gzip_cmd, "unpaired_")
            up_read_count, chnk_num, up_sp = t
            val = [prev_line[0], prev_line[1], prev_line[2]]
            up_sp.stdin.write(up_outstring.format(*val))
            up_read_count += 1
            tot_read_count += 1
            prev_line = [bam_line[0]] + bam_line[9:11] + sam_flag_list[0:3]
    else:
        bam_infile.close()

        if len(prev_line):
            t = check_read_count(up_read_count, chnk_num, outdir,
                                 chnk_size, up_sp, gzip_cmd, "unpaired_")
            up_read_count, chnk_num, up_sp = t
            val = [prev_line[0], prev_line[1], prev_line[2]]
            up_sp.stdin.write(up_outstring.format(*val))
            up_read_count += 1
            tot_read_count += 1

        write_interval_file(outdir)
        print tot_read_count


def decompose_flag(flag, bits):
    """Identifies the bits that are set in the SAM flag"""
    flag_list = [bool(int(flag) & x) for x in bits]
    return flag_list


def check_read_count(rc, chnk_num, chnk_dir, chnk_size, sp, cmd, read_type):
    """Checks if read count limit reached, spawns new gzip process."""
    if rc >= chnk_size:
        chnk_num += 1
        sp = spawn_gzip(chnk_num, chnk_dir, cmd, read_type)
        rc = 0
    return [rc, chnk_num, sp]


def revcomp(seq):
    """Reverse complements a sequence"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rev = "".join(complement.get(base, base) for base in reversed(seq))
    return rev


def spawn_gzip(chnk_num, chnk_dir, gzip_cmd, read_type):
    """Returns gzip process filehandle."""
    s = os.path.join(chnk_dir, str(read_type) + "chunk" + str(chnk_num) + '.fastq.gz')
    p = subprocess.Popen(gzip_cmd + s, bufsize=-1, shell=True,
                         stdin=subprocess.PIPE)
    return p


def write_interval_file(directory):
    """Creates an interval file in the given directory"""
    interval_filename = os.path.join(directory, 'interval.txt')
    outstring = [ str(x) + '\n' for x in os.listdir(directory)
                 if os.path.isfile(os.path.join(directory, x))]
    interval_file = open(interval_filename, 'w')
    interval_file.write(string.join(outstring, ''))
    interval_file.close()


def parse_args():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--num_reads', '-n', type=int, default=5000000000,
                        help='Specify number of reads per FASTQ file.')
    parser.add_argument('bam', nargs=1, type=argparse.FileType('r'),
                        help='Specify SAM file')
    parser.add_argument('outdir', default='.', nargs=1)
    return parser.parse_args()


if __name__ == '__main__':
    main()
