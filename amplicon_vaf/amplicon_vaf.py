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

# Add augment_maf to sys.path
sys.path.insert(1, os.path.join(sys.path[0], "..", "augment_maf"))
from augment_maf import bamUtils
from augment_maf import mafUtils

def main():
    #get options
    parser = argparse.ArgumentParser(description='obtain input files and mutation details, if available')
    parser.add_argument('-b', '--bam', nargs='+', type=str, dest='bam_filenames', help='bam file for sample provided')
    parser.add_argument('-m', '--maf', dest='maf', action='append', default = [], help='MAF file(s) containing details for mutations of interest')
    parser.add_argument('-r', '--ref', dest='genome', help='indexed fasta file alinged to')
    parser.add_argument('-o', '--output', dest='outfile', help='output file')
    parser.add_argument('-M', '--mode', default='pileup', choices=bamUtils.MODES.keys())
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

    count_keys = ["Reference_Count", "Alternate_Count", "VAF", "Sample_ID", "Variant_ID"]
    fields = readers[0].fieldnames
    for key in count_keys:
        if not key in fields:
            fields.append(key)

    if args.outfile:
        outfile = open(args.outfile, "w")
    else:
        outfile = sys.stdout
    writer = csv.DictWriter(outfile, delimiter="\t", lineterminator="\n", fieldnames=fields)

    bams = [pysam.AlignmentFile(bam) for bam in args.bam_filenames]
    # Take bam filename up until first '.' as sample id
    sample_ids = [os.path.basename(bam).split(os.path.extsep)[0] for bam in args.bam_filenames]

    ref_key = "Reference_Count"
    alt_key = "Alternate_Count"
    vaf_key = "VAF"
    sample_id_key = "Sample_ID"
    variant_id_key = "Variant_ID"
    
    writer.writeheader()
    for (sample_id, bam) in zip(sample_ids, bams):
        for row in reader:
            chrom = chr_prefix + row["Chromosome"]
            pos = int(row["Start_Position"])
            ref = row["Reference_Allele"]
            alt = mafUtils.get_nref_allele(row)
            
            if mafUtils.is_snv(row):
                counts = bamUtils.count_bases(bam, reffile, chrom, pos)
            else:
                counts = bamUtils.count_indels(bam, reffile, chrom, pos, ref, alt, args.mode)
            row[ref_key] = counts[ref]
            row[alt_key] = counts[alt]

            try:
                vaf = counts[alt] / (counts[alt] + counts[ref])
            except ZeroDivisionError:
                vaf = float('nan')

            if (math.isnan(vaf)):
                row[vaf_key] = "NaN"
            else:
                row[vaf_key] = round(vaf, 6)
                

            row[sample_id_key] = sample_id

            try:
                variant_id = row[variant_id_key]
            except KeyError:
                variant_id = str(chrom) + ":" + str(pos) + ref + ">" + alt

            row[variant_id_key] = variant_id
            
            writer.writerow(row)

if __name__ == "__main__":
    main()
