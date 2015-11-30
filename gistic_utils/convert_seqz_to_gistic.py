"""
convert_seqz_to_gistic.py
==========================
Description
Converts Sequenza segmentation file into a GISTIC segmentation file format.
Adds number of markers per sequenza segment with windows of specified size
(default: 100,000bp). Outputs to STDOUT.

Assumes the first line in the Sequenza segmentation file is a header line.

Dependencies:
-None

Inputs
-Sequenza segmentation file
-Sample ID for Sequenza segmentation file

Outputs
-Prints GISTIC formatted segmentation file to STDOUT

Known Issues
-None
"""

import argparse
import math
import sys


def main():
    args = parse_args()
    window_size = args.window_size
    seqz_seg = args.seqz_seg
    sample_id = args.sample_id
    header_found = False
    #print "Sample\tChromosome\tStart Position\tEnd Position\tNum markers\tSeg.CN"
    for line in seqz_seg:
        if not header_found:
            header_found = True
            continue
        fields = line.split()
        chrm = fields[0]
        if chrm.startswith('"') and chrm.endswith('"'):
            chrm = chrm[1:-1]
        start = int(fields[1])
        end = int(fields[2])
        cnt = int(fields[9])
        s = int(start/window_size) * window_size
        seg_total = 0
        while s < end:
            seg_total += 1
            s += window_size
        # convert CNt to log2
        if cnt > 0:
            cnt = math.log(cnt/float(2), 2)
        elif cnt == 0:
            cnt = -1
        else:
            sys.stderr.write("CNt {} is less than 0.".format(cnt))
        print "{}\t{}\t{}\t{}\t{}\t{}".format(sample_id, chrm,
                                              start, end, seg_total, cnt)
    return


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqz_seg', type=argparse.FileType('r'),
                        help='Sequenza segment file to convert.')
    parser.add_argument('--window_size', '-w', type=int, default=100000,
                        help='Size of window marker.')
    parser.add_argument('--sample_id', '-s', required=True,
                        help='Specify sampleID for GISTIC seg file.')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
