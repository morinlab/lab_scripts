#!/usr/bin/env python

"""Variant calling from Access Array data to extract the VAF of known mutations."""
from __future__ import division
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
    #get options
    parser = argparse.ArgumentParser(description='obtain input files and mutation details, if available')
    parser.add_argument('-b', '--bam', nargs='+', type=str, dest='bam_filenames', help='bam file for sample provided')
    parser.add_argument('-m', '--maf', dest='maf', action='append', default = [], help='MAF file(s) containing details for mutations of interest')
    parser.add_argument('-r', '--ref', dest='genome', help='indexed fasta file alinged to')
    parser.add_argument('-o', '--output', dest='outfile', help='output file')
    args = parser.parse_args()
    
    try:
        reffile = pysam.Fastafile(args.genome)
    except AttributeError:
        print "Error finding reference file:" + args.genome
	sys.exit(-1)

    if reffile.references[0].startswith("chr"):
        chr_prefix = "chr"
    else:
        chr_prefix = ""
    
    readers = []
    for f in args.maf:
        lines = [x for x in open(f) if not x.startswith("#")]
        readers.append(csv.DictReader(lines, delimiter="\t"))
    reader = list(*readers)

    count_keys = ["ref_count", "alt_count", "vaf", "sample_id"]
    fields = readers[0].fieldnames
    for key in count_keys:
        if not key in fields:
            fields.append(key)

    if args.outfile:
        outfile = open(args.outfile, "w")
    else:
        outfile = sys.stdout
    writer = csv.DictWriter(outfile, delimiter="\t", fieldnames=fields)
    # writer.writeheader()

    bams = [pysam.AlignmentFile(bam) for bam in args.bam_filenames]
    # Take bam filename up until first '.' as sample id
    sample_ids = [os.path.basename(bam).split(os.path.extsep)[0] for bam in args.bam_filenames]

    ref_key = "ref_count"
    alt_key = "alt_count"
    vaf_key = "vaf"
    sample_id_key = "sample_id"

    for (sample_id, bam) in zip(sample_ids, bams):
        for row in reader:
            chrom = chr_prefix + row["Chromosome"]
            pos = int(row["Start_Position"])
            ref = row["Reference_Allele"]
            alt = mafUtils.get_nref_allele(row)

            if mafUtils.is_snv(row):
                counts = bamUtils.count_bases(bam, reffile, chrom, pos)
            else:
                counts = bamUtils.count_indels(bam, reffile, chrom, pos, ref, alt, "hybrid")
     
            row[ref_key] = counts[ref]
            row[alt_key] = counts[alt]

            vaf = counts[alt] / (counts[alt] + counts[ref])
            row[vaf_key] = round(vaf, 6)
            row[sample_id_key] = sample_id
            writer.writerow(row)

if __name__ == "__main__":
    main()
