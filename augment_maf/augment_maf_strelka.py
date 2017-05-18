#!/usr/bin/env python2

import argparse
import warnings
import csv

import vcf as pyvcf

import mafUtils
import strelkaUtils


def main(args):
    # Create MAF reader
    lines = [x for x in open(args.maf) if not x.startswith("#")]
    maf_reader = csv.DictReader(lines, delimiter="\t")

    # Confirm existence of count keys and update accordingly
    count_keys = ["n_ref_count", "n_alt_count",
                  "t_ref_count", "t_alt_count",
                  "t_depth", "n_depth"]
    fields = maf_reader.fieldnames
    for key in count_keys:
        if key not in fields:
            fields.append(key)
    fields.append("strelka_qual")

    # Create MAF writer
    writer = csv.DictWriter(args.output, delimiter="\t", fieldnames=fields)
    writer.writeheader()

    # Create index of Strelka allele counts
    vcf_reader = pyvcf.Reader(filename=args.vcf)
    count_idx = strelkaUtils.generate_count_index(vcf_reader)

    # Iterate over rows in the MAF file
    for row in maf_reader:
        # Extract row features
        chrom = row["Chromosome"]
        pos = int(row["Start_Position"])
        ref = row["Reference_Allele"]
        idx_key = (chrom, pos, ref)
        old_counts = mafUtils.get_allele_counts(row)
        # Reset allele counts, if applicable
        if sum(old_counts) > 0 and args.replace:
            old_counts = (0, 0, 0, 0)
        row.update(dict(zip(count_keys, old_counts)))
        # Extract allele counts from Strelka index
        if idx_key in count_idx:
            t_ref_count, t_alt_count, n_ref_count, n_alt_count, strelka_qual = count_idx[idx_key]
        else:
            warnings.warn("Variant {} not found in Strelka VCF file".format(idx_key))
            continue
        # Update allele counts
        row["t_ref_count"] += t_ref_count
        row["t_alt_count"] += t_alt_count
        row["t_depth"] = row["t_ref_count"] + row["t_alt_count"]
        row["n_ref_count"] += n_ref_count
        row["n_alt_count"] += n_alt_count
        row["n_depth"] = row["n_ref_count"] + row["n_alt_count"]
        # Add Strelka quality information
        row["strelka_qual"] = strelka_qual
        # Output updated row
        writer.writerow(row)


def parse_args():
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("maf", help="MAF file")
    parser.add_argument("vcf", help="Associated Strelka VCF file")
    parser.add_argument("--replace", "-r", action="store_true",
                        help="Replace instead of adding to the existing values")
    parser.add_argument("--output", "-o", type=argparse.FileType("w"))
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    main(args)
