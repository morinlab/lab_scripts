#!/usr/bin/env python

"""
calc_capture_stats
==================
This script calculates useful statistics for capture-based experiments.
Notably, it calculates:
- Average genome coverage
- Average target coverage
- Percent on-target reads
- Percent fold enrichment

Inputs:
- BAM file (sorted, rmduped and indexed; only one read length)
- A set of intervals (BED file)

Output:
- Simple TSV file with stats
"""

import argparse
import cancer_api
import pysam
from collections import OrderedDict


def main():

    # Argument parsing
    parser = argparse.ArgumentParser(description="Calculate useful stats for capture-based "
                                     "experiments.")
    parser.add_argument("bam_file", help="Input BAM file (sorted, rmduped and indexed). "
                        "Assumes that the BAM file contains reads of one length.")
    parser.add_argument("bed_file", help="Input BED intervals file.")
    parser.add_argument("--output", help="Output TSV stats file.")
    args = parser.parse_args()

    # Initialize variables
    inbam = pysam.AlignmentFile(args.bam_file, "rb")
    inbed = cancer_api.files.BedFile.open(args.bed_file)
    output = args.output if args.output else inbam.filename + ".capture_stats.tsv"
    stats = OrderedDict()

    # Calculate overall coverage
    read_length = inbam.head(1).next().query_length
    genome_num_mapped = inbam.mapped
    genome_length = sum(inbam.lengths)
    genome_cov = read_length * genome_num_mapped / float(genome_length)
    stats["Genome_Coverage"] = round(genome_cov, 2)

    # Calculate on-target coverage
    target_length = 0
    target_num_mapped = 0
    for interval in inbed:
        target_length += interval.length
        reads_iterator = inbam.fetch(interval.chrom, interval.start_pos, interval.end_pos)
        target_num_mapped += len([1 for r in reads_iterator])
    target_cov = read_length * target_num_mapped / float(target_length)
    stats["Target_Coverage"] = round(target_cov, 2)

    # Calculate percent on-target
    percent_on_target = target_num_mapped / float(genome_num_mapped) * 100
    stats["Percent_On_Target"] = round(percent_on_target, 2)

    # Calculate percent fold enrichment
    percent_fold_enrichment = target_cov / genome_cov * 100
    stats["Percent_Fold_Enrichment"] = round(percent_fold_enrichment, 2)

    # Write out stats to file
    with open(output, "w") as outfile:
        for k, v in stats.items():
            stats[k] = str(v)
        outfile.write("\t".join(stats.keys()) + "\n")
        outfile.write("\t".join(stats.values()) + "\n")

    # Cleanup
    inbam.close()

if __name__ == '__main__':
    main()
