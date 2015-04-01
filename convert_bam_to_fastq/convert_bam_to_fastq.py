#!/usr/bin/env python

"""
convert_bam_to_fastq.py
=======================
This script converts a BAM file into a FASTQ file.

Inputs:
- BAM file

Outputs:
- FASTQ file

Requirements
------------
- cancer_api >= v0.1.6
- pysam >= v0.8.1
- The input BAM file should contain unpaired reads.
  This script was created for the realignment pipeline,
  where unpaired reads need to be converted into a FASTQ
  file before being realigned. Accordingly, we expect read
  names to appear uniquely in the BAM file. If not, only
  one of the two will be outputted.

Known Issues
------------
- None
"""

import argparse
import logging
import cancer_api
import pysam

__version__ = "v1.0.0"


def main():

    # ========================================================================================== #
    # Argument parsing
    # ========================================================================================== #

    parser = argparse.ArgumentParser(description="Convert a BAM file into FASTQ file(s)")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_fastq", help="Output FASTQ file")
    args = parser.parse_args()

    # ========================================================================================== #
    # Set up
    # ========================================================================================== #

    # Setup logging
    cancer_api.utils.setup_logging()

    # Instantiate read buffer
    added_reads = {}

    # Open various files
    inbam = pysam.AlignmentFile(args.input_bam, "rb")
    outfastq = cancer_api.FastqFile.new(args.output_fastq, buffersize=5000000)

    # ========================================================================================== #
    # Convert BAM file into FASTQ file(s)
    # ========================================================================================== #

    # Iterate over reads in BAM file and partition into R1 and R2
    logging.info("Iterating over reads in BAM file and adding to FASTQ file...")
    for read in inbam:
        if read.query_name in added_reads:
            logging.warn("Skipping read because read name has already come up ({})...".format(
                         read.query_name))
            continue
        rawread = cancer_api.misc.RawRead(
            read.query_name, read.query_sequence, "+",
            "".join([str(unichr(x + 33)) for x in read.query_qualities]))
        outfastq.add_obj(rawread)
        added_reads[read.query_name] = 1
    # Ensure that any remaining reads are written to disk
    outfastq.write()


if __name__ == '__main__':
    main()
