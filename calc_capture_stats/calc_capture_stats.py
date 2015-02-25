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
- A set of excluded intervals (optional; BED file)

Output:
- Simple TSV file with stats
"""

import argparse
import re
import cancer_api
import pysam
from collections import OrderedDict

__version__ = "1.1.0"


def main():

    # Argument parsing
    parser = argparse.ArgumentParser(description="Calculate useful stats for "
                                     "capture-based experiments.")
    parser.add_argument("bam_file", help="Input BAM file (sorted, rmduped and indexed). "
                        "Assumes that the BAM file contains reads of one length.")
    parser.add_argument("bed_file", help="BED file listing target intervals.")
    parser.add_argument("--output", help="Output TSV stats file.")
    parser.add_argument("--exclude", help="BED file listing regions to be excluded. "
                        "Useful for pooled libraries sharing the same multiplex index.")
    args = parser.parse_args()

    # Initialize variables
    inbam = pysam.AlignmentFile(args.bam_file, "rb")
    inbed = cancer_api.files.BedFile.open(args.bed_file)
    excl_bed = cancer_api.files.BedFile.open(args.exclude) if args.exclude else None
    output = args.output if args.output else inbam.filename + ".capture_stats.tsv"
    stats = OrderedDict()

    # Define helper functions
    def count_interval_reads(bam_file, interval, viewed_reads={}):
        """Counts the number of reads in a given interval in a BAM file.
        Adds viewed reads to dictionary to avoid double-counting of reads.
        """
        reads_iterator = inbam.fetch(interval.chrom, interval.start_pos, interval.end_pos)
        read_qnames = [r.query_name for r in reads_iterator
                       if r.query_name not in viewed_reads and not r.is_duplicate]
        num_reads = len(read_qnames)
        for qname in read_qnames:
            viewed_reads[qname] = 1
        return num_reads

    # Calculate overall coverage
    read_length = inbam.head(1).next().query_length
    # Using flagstat in order to not count duplicate reads
    total_index = 0
    dups_index = 3
    flagstat_regex = r"(\d+).*"
    flagstat = pysam.flagstat(inbam.filename)
    num_total = int(re.match(flagstat_regex, flagstat[total_index]).group(1))
    num_dups = int(re.match(flagstat_regex, flagstat[dups_index]).group(1))
    genome_num_mapped = num_total - num_dups
    # Calculate genome length
    genome_length = sum(inbam.lengths)
    if excl_bed:
        # Correct for excluded regions
        excl_length = 0
        excl_num_mapped = 0
        for interval in excl_bed:
            excl_length += interval.length
            excl_num_mapped += count_interval_reads(inbam, interval)
        genome_num_mapped -= excl_num_mapped
        genome_length -= excl_length
    genome_cov = read_length * genome_num_mapped / float(genome_length)
    stats["Genome_Coverage"] = round(genome_cov, 3)

    # Calculate on-target coverage
    target_length = 0
    target_num_mapped = 0
    viewed_reads = {}
    for interval in inbed:
        target_length += interval.length
        target_num_mapped += count_interval_reads(inbam, interval, viewed_reads)
    target_cov = read_length * target_num_mapped / float(target_length)
    stats["Target_Coverage"] = round(target_cov, 3)

    # Calculate percent on-target
    percent_on_target = target_num_mapped / float(genome_num_mapped) * 100
    stats["Percent_Reads_On_Target"] = round(percent_on_target, 3)

    # Calculate percent fold enrichment
    percent_fold_enrichment = target_cov / genome_cov * 100
    stats["Percent_Fold_Enrichment"] = round(percent_fold_enrichment, 3)

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
