#!/usr/bin/env python

# This script augments the reference and normal counts in a MAF file with the
# counts found in one or more BAM files. At least one MAF file must be
# specified. If more than one MAF file is specified, the output will be a
# pooled MAF file with the unique entries from all the input files.

import pysam
import csv
import argparse
import mafUtils
import bamUtils
import warnings
import sys
import itertools

if __name__ == "__main__":
    desc = "Add normal and/or tumour allele counts to a MAF file"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--replace", "-r", action="store_true", default=False)
    parser.add_argument("--normal-bam", "-n", action="append", default=[])
    parser.add_argument("--tumour-bam", "-t", action="append", default=[])
    parser.add_argument("--maf", "-m", action="append", default=[])
    parser.add_argument("reference")
    parser.add_argument("outfile")
    args = parser.parse_args()

    if (args.maf == []):
        sys.exit("You must specify at least one MAF file")

    reffile = pysam.Fastafile(args.reference)
    normal_sams = [pysam.Samfile(bam) for bam in args.normal_bam]
    tumour_sams = [pysam.Samfile(bam) for bam in args.tumour_bam]
    samfiles = {"normal": normal_sams, "tumour": tumour_sams}
    
    readers = []
    for f in args.maf:
        lines = [x for x in open(f) if not x.startswith("#")]
        readers.append(csv.DictReader(lines, delimiter="\t"))
    reader = itertools.chain(*readers)

    count_keys = ["n_ref_count", "n_alt_count", "t_ref_count", "t_alt_count"]
    fields = readers[0].fieldnames
    for key in count_keys:
        if not key in fields:
            fields.append(key)

    outfile = open(args.outfile, "w")
    writer = csv.DictWriter(outfile, delimiter="\t", fieldnames=fields)
    writer.writeheader()

    done = []
    for row in reader:
        chrom = row["Chromosome"]
        pos = int(row["Start_Position"])
        ref = row["Reference_Allele"]
        alt = mafUtils.get_nref_allele(row)
        if (chrom, pos, ref, alt) in done:
            continue

        old_counts = mafUtils.get_allele_counts(row)
        if sum(old_counts) > 0:
            if args.replace:
                warnings.warn("Replacing existing allele counts")
                old_counts = (0, 0, 0, 0)
            else:
                warnings.warn("Adding to existing allele counts")
        
        row.update(dict(zip(count_keys, old_counts)))

        for sample, sams in samfiles.items():
            ref_key = "{}_ref_count".format(sample[0])
            alt_key = "{}_alt_count".format(sample[0])
            for samfile in sams:
                if mafUtils.is_snv(row):
                    counts = bamUtils.count_bases(samfile, reffile, chrom, pos)
                else:
                    counts = bamUtils.count_indels(samfile, reffile, chrom, pos, ref, alt)
                row[ref_key] += counts[ref]
                row[alt_key] += counts[alt]

        writer.writerow(row)
        done.append((chrom, pos, ref, alt))
