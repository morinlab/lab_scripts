#!/usr/bin/env python

"""
calc_stelka_vaf
===============
This script calculates the VAF of variants called
by Strelka.
"""

from __future__ import division

import sys
import argparse
import warnings

import vcf as pyvcf


def main():

    # Parse command-line arguments
    args = parse_args()

    # Iterate over VCF files
    for vcf in args.vcfs:
        with open(vcf) as ovcf:
            vcf_reader = pyvcf.Reader(ovcf)
            for record in vcf_reader:
                if record.is_snp:
                    ref_count, alt_count = calc_snv_counts(record)
                elif record.is_indel:
                    ref_count, alt_count = calc_indel_counts(record)
                else:
                    warnings.warn("Encountered variant that isn't a SNV or an indel.")
                args.output.write(output_vaf(record, ref_count, alt_count))


def parse_args():
    """Parse command-line arguments"""
    # Initialize parser
    parser = argparse.ArgumentParser()
    # Positional arguments
    parser.add_argument("vcfs", metavar="vcf_file", nargs="+", help="Strelka VCF file(s)")
    # Optional arguments
    parser.add_argument("--output", "-o", type=argparse.FileType("r"), default=sys.stdout, help="Output file")
    # Parse arguments and return
    args = parser.parse_args()
    return args


def calc_snv_counts(record):
    # Create keys to extract read counts
    ref_key = "{}U".format(record.REF)
    alt_key = "{}U".format(record.ALT[0])
    ref_supp = record.samples[1][ref_key]
    alt_supp = record.samples[1][alt_key]
    # Only consider tier 1 reads
    ref_count = int(ref_supp[0])
    alt_count = int(alt_supp[0])
    return (ref_count, alt_count)


def calc_indel_counts(record):
    """Calculate the number of reads supporting the
    reference and alternate alleles based on what's
    in the Strelka VCF file.
    """
    # Get reads supporting "alternate" (ref) and indel alleles
    tar = record.samples[1]["TAR"]
    tir = record.samples[1]["TIR"]
    # Only consider tier 1 reads
    ref_count = int(tar[0])
    alt_count = int(tir[0])
    return (ref_count, alt_count)


def output_vaf(record, ref_count, alt_count):
    """Output string for output file for a given
    variant.
    """
    cols = []
    cols.extend([record.CHROM, record.POS, record.REF, record.ALT[0]])
    cols.extend([ref_count, alt_count])
    vaf = round(alt_count / (ref_count + alt_count), 2)
    cols.append(vaf)
    output = "\t".join(map(str, cols)) + "\n"
    return output


if __name__ == '__main__':
    main()
