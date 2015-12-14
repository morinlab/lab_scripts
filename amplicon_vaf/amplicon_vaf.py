#!/usr/bin/env python

"""Variant calling from Access Array data to extract the VAF of known mutations."""
import pysam
import csv
import itertools
import math
import sys
import argparse
import os.path
import warnings

from augment_maf import bamUtils
from augment_maf import mafUtils

def main():
    print "DEBUG: Using pysam version: " + pysam.__version__ 
    #set some global variables

    #get options
    parser = argparse.ArgumentParser(description='obtain input files and mutation details, if available')
    parser.add_argument('-s', '--sample', dest='sample', help='name for sample being provided')
    parser.add_argument('-b', '--bam', nargs='+', type=str, dest='bam_filenames', help='bam file for sample provided')
    parser.add_argument('-m', '--maf', dest='maf', action='append', default = [], help='MAF file(s) containing details for mutations of interest')
    parser.add_argument('-r', '--region_file', dest='regions', help='file containing coordinates of one or more regions targeted in sample(s)')
    parser.add_argument('-g', '--genome_fasta', dest='genome', help='indexed fasta file alinged to')
    parser.add_argument('-o', '--output', dest='outfile', help='output file')
    args = parser.parse_args()
    
    try:
        reffile = pysam.Fastafile(args.genome)
    except AttributeError:
        print "using default genome fasta %s" % genome_fasta

    if reffile.references[0].startswith("chr"):
        chr_prefix = "chr"
    else:
        chr_prefix = ""
    
    readers = []
    for f in args.maf:
        lines = [x for x in open(f) if not x.startswith("#")]
        readers.append(csv.DictReader(lines, delimiter="\t"))
    reader = itertools.chain(*readers)

    count_keys = ["ref_count", "alt_count"]
    fields = readers[0].fieldnames
    for key in count_keys:
        if not key in fields:
            fields.append(key)

    if args.outfile:
        outfile = open(args.outfile, "w")
    else:
        outfile = sys.stdout
    writer = csv.DictWriter(outfile, delimiter="\t", fieldnames=fields)
    writer.writeheader()

    bams = [pysam.AlignmentFile(bam) for bam in args.bam_filenames]
    bamfiles = {"sample": bams}
    
    done = []
    for row in reader:
        chrom = chr_prefix + row["Chromosome"]
        pos = int(row["Start_Position"])
        ref = row["Reference_Allele"]
        alt = mafUtils.get_nref_allele(row)
        if (chrom, pos, ref, alt) in done:
            continue

        old_counts = mafUtils.get_allele_counts(row)

        if sum(old_counts) > 0:
            if args.replace:
                warnings.warn("Replacing existing allele counts")
                old_counts = (0, 0)
            else:
                warnings.warn("Adding to existing allele counts")

        row.update(dict(zip(count_keys, old_counts)))

        for bam in bams:
            ref_key = "ref_count"
            alt_key = "alt_count"
            if mafUtils.is_snv(row):
                counts = bamUtils.count_bases(bam, reffile, chrom, pos)
            else:
                counts = bamUtils.count_indels(bam, reffile, chrom, pos, ref, alt, "hybrid")
            row[ref_key] += counts[ref]
            row[alt_key] += counts[alt]

        writer.writerow(row)
        done.append((chrom, pos, ref, alt))

if __name__ == "__main__":
    main()
