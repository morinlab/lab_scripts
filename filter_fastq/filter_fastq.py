#!/usr/bin/env python

"""
filter_fastq.py
==============
This script filters FASTQ files according to specified
regular expressions for read sequences.

Inputs:
- FASTQ file(s) (single- or paired-end)
- Regular expression(s) for filtering

Outputs:
- Filtered FASTQ file(s)

Requirements
------------
- cancer_api >= 0.2.0

Known Issues
------------
- In the current implementation, if you specify both an include and
    an exclude regular expression, the exclude filter takes precedence.
"""

import argparse
import os
import logging
import re
import cancer_api

__version__ = "v1.0.0"


def main():

    # ========================================================================================== #
    # Argument parsing
    # ========================================================================================== #

    parser = argparse.ArgumentParser(description="Filter FASTQ files based on an include and/or "
                                     "an exclude regular expression(s).")
    parser.add_argument("fastq", nargs="+", help="FASTQ file(s) (single- or paired-end)")
    parser.add_argument("--include_regex", "-i", help="Regular expression for included reads")
    parser.add_argument("--exclude_regex", "-e", help="Regular expression for excluded reads")
    parser.add_argument("--output_dir", default=".", help="Output directory")
    parser.add_argument("--num_buffer", "-b", type=int, default=2500000,
                        help="Number of reads kept in buffer before flushing to disk")
    args = parser.parse_args()

    # ========================================================================================== #
    # Set up
    # ========================================================================================== #

    # Setup logging
    cancer_api.utils.setup_logging()

    # Check if one or two FASTQ files
    if len(args.fastq) == 1:
        fastq1_filepath = args.fastq[0]
        fastq2_filepath = None
    elif len(args.fastq) == 2:
        fastq1_filepath = args.fastq[0]
        fastq2_filepath = args.fastq[1]
    else:
        raise ValueError("Did not receive one or two FASTQ files.")

    # Check if output_dir exists; if not, create it
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(args.output_dir)

    # Setting up files
    logging.info("Opening FASTQ files...")
    infastq1 = cancer_api.files.FastqFile.open(fastq1_filepath)
    root, ext = infastq1.split_filename()
    outfastq1_filepath = root + ".filtered." + ext
    outfastq1 = cancer_api.files.FastqFile.new(
        os.path.join(output_dir, outfastq1_filepath), buffersize=args.num_buffer)
    infastq2 = None
    if fastq2_filepath:
        infastq2 = cancer_api.files.FastqFile.open(fastq2_filepath)
        root, ext = infastq2.split_filename()
        outfastq2_filepath = root + ".filtered." + ext
        outfastq2 = cancer_api.files.FastqFile.new(
            os.path.join(output_dir, outfastq2_filepath), buffersize=args.num_buffer)

    # More variables
    iregex = re.compile(args.include_regex) if args.include_regex else args.include_regex
    eregex = re.compile(args.exclude_regex) if args.exclude_regex else args.exclude_regex

    # Determine default filtering status
    if iregex and not eregex:
        default_is_filtered = True
    elif not iregex and eregex:
        default_is_filtered = False
    elif iregex and eregex:
        default_is_filtered = True
    elif not iregex and not eregex:
        raise ValueError("No regular expression given for filtering.")

    # ========================================================================================== #
    # Filter FASTQ files
    # ========================================================================================== #

    # Iterate over reads
    logging.info("Filtering reads in FASTQ files...")
    if infastq2:
        infastq2_iter = infastq2.__iter__()
    for read1 in infastq1:
        seq1 = read1.seq
        if infastq2:
            read2 = next(infastq2_iter)
            seq2 = read2.seq
        is_filtered1 = default_is_filtered
        # Apply include regex
        if iregex and iregex.search(seq1):
            is_filtered1 = False
        # Apply exclude regex
        if eregex and eregex.search(seq1):
            is_filtered1 = True
        # Do the same for FASTQ #2, if applicable
        is_filtered2 = default_is_filtered
        if infastq2:
            # Apply include regex
            if iregex and iregex.search(seq2):
                is_filtered2 = False
            # Apply exclude regex
            if eregex and eregex.search(seq2):
                is_filtered2 = True
        # Only write out non-filtered reads
        # Check if one of the is_filtered variables has changed.
        # If so, use the changed value to determine whether to filter
        final_is_filtered = default_is_filtered
        if is_filtered1 is not default_is_filtered:
            final_is_filtered = is_filtered1
        elif is_filtered2 is not default_is_filtered:
            final_is_filtered = is_filtered2
        if not final_is_filtered:
            outfastq1.add_obj(read1)
            if infastq2:
                outfastq2.add_obj(read2)

    # Write out any remaining reads
    outfastq1.write()
    if infastq2:
        outfastq2.write()


if __name__ == '__main__':
    main()
