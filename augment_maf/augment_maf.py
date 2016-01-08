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
import logging

if __name__ == "__main__":
    desc = "Add normal and/or tumour allele counts to a MAF file"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--replace", "-r", action="store_true", default=False)
    parser.add_argument("--normal-bam", "-n", action="append", default=[])
    parser.add_argument("--tumour-bam", "-t", action="append", default=[])
    parser.add_argument("--maf", "-m", action="append", default=[])
    parser.add_argument("--mode", default="hybrid", choices=bamUtils.MODES.keys())
    parser.add_argument("--log_lvl", type=str.upper, default="INFO")
    parser.add_argument("--log_file", type=argparse.FileType("w"), default=sys.stderr)
    parser.add_argument("reference")
    parser.add_argument("outfile")
    args = parser.parse_args()

    # Setup logging
    log_lvl = getattr(logging, args.log_lvl)
    log_fmt = '%(asctime)s - %(levelname)s (%(module)s.%(funcName)s):  %(message)s'
    date_fmt = '%Y/%m/%d %H:%M:%S'  # 2010/12/12 13:46:36
    logging.basicConfig(stream=args.log_file, level=log_lvl, format=log_fmt, datefmt=date_fmt)

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
        if key not in fields:
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
                    counts = bamUtils.count_indels(samfile, reffile, chrom, pos, ref, alt, args.mode)
                row[ref_key] += counts[ref]
                row[alt_key] += counts[alt]
            # Update depth columns as well
            depth_key = "{}_depth".format(sample[0])
            row[depth_key] = row[ref_key] + row[alt_key]

        writer.writerow(row)
        done.append((chrom, pos, ref, alt))
    if args.log_file:
        args.log_file.close()
