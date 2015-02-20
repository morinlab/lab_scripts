#!/usr/bin/env python

"""
filter_vcf_by_reads.py
======================
This script filters the variants listed in a VCF file
according to a minimum number of supporting reads.
This is done by inspecting the matched BAM file.

Known Issues
------------
- We need to settle on a sensible default for the
    minimum number of reads supporting a variant.
"""

import argparse
import pysam


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input_vcf', help='Input VCF file')
    parser.add_argument('input_bam', help='BAM file matched to input VCF file')
    parser.add_argument('--output_vcf', help='Filtered VCF file (default: *.filtered.vcf)')
    parser.add_argument('--read_cutoff', default=10,
                        help='Minimum number of unique read IDs required to keep a variant')
    args = parser.parse_args()

    print args

    # Open VCFs for reading and writing
    fh_in_vcf = open(args.in_vcf, 'r')
    fh_out_vcf = open(args.out_vcf, 'w')

    # Open the BAM file to be used for extracting read data
    bamfile = pysam.Samfile(args.bam_file, "rb")

    # Keep track of how many variants we have removed
    rm_counter = 0

    # Now process each line in the input VCF, looking for the amount of read support for each variant
    for line in fh_in_vcf:
        # If this is a header line, print to output VCF immediately
        if line[0] == "#":
            fh_out_vcf.write(line)
            continue

        # Otherwise, extract info about this variant
        line_split = line.split("\t")
        chrom = line_split[0]
        pos = int(line_split[1])
        ref = line_split[3]
        alt = line_split[4]

        # If this is an indel, skip it
        if len(ref) > 1 or len(alt) > 1:
            fh_out_vcf.write(line)
            continue

        # Create a pileup for reads surrounding the variant
        var_reads = []
        for pcol in bamfile.pileup(reference=chrom, start=pos-1, end=pos+1):
            # Look at the pileup data for the specific position we are interested in
            genome_coord = pcol.pos + 1
            if genome_coord != pos:
                continue

            for pread in pcol.pileups:
                if pread.alignment.seq[pread.qpos] == alt:
                    if pread.alignment.qname not in var_reads:
                        var_reads.append(pread.alignment.qname)

        # Check if we have found the minimum number of reads supporting the variant
        if len(var_reads) >= int(args.read_cutoff):
            fh_out_vcf.write(line)
        else:
            print "Removed %s:%s -> %s reads" % (chrom, str(pos), str(len(var_reads)))
            rm_counter = rm_counter + 1

    print "\nTOTAL REMOVED = %s" % str(rm_counter)


if __name__ == '__main__':
    main()
