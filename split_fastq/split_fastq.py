#!/usr/bin/env python

"""
split_fastq.py
==============
This script splits FASTQ files into smaller ones.
This is meant to allow for parallel alignment of
large NGS datasets (e.g., genome data).

Inputs:
- FASTQ file(s) (single- or paired-end)

Outputs:
- Output directory into which all smaller FASTQ files
    will be stored

Known Issues
------------
- None
"""

import argparse
import os
import logging
import cancer_api

__version__ = "v1.1.3"

# MIN_SPILLOVER indicates the minimum fraction (between 0 and 1) of num_reads
# that is required to create a new chunk. This is mostly meant to prevent the
# creation of small FASTQ files with only very few reads.
MIN_SPILLOVER = 0.25


def main():

    # ========================================================================================== #
    # Argument parsing
    # ========================================================================================== #

    parser = argparse.ArgumentParser(description="Split FASTQ files into smaller ones.")
    parser.add_argument("fastq", nargs="+", help="FASTQ file(s) (single- or paired-end)")
    parser.add_argument("--num_reads", "-n", type=int, default=75000000,
                        help="Number of reads per output FASTQ file")
    parser.add_argument("--num_buffer", "-b", type=int, default=5000000,
                        help="Number of reads kept in buffer before flushing to disk")
    parser.add_argument("--output_dir", default=".", help="Output directory")
    parser.add_argument("--interval_file", help="Output intervals (for use in Pipeline Factory)")
    parser.add_argument("--no_compression", action="store_true",
                        help="Disables gzip compression of output FASTQ files")
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

    # ========================================================================================== #
    # Split FASTQ files
    # ========================================================================================== #

    # Define helper classes and functions
    class Chunk(object):
        """Convenience class for managing chunks"""
        def __init__(self, infastq, output_dir, no_compression):
            """Initialize first chunk"""
            self.infastq = infastq
            self.filename = infastq.filename
            self.output_prefix = os.path.join(output_dir, self.filename)
            self.counter = 0
            self.outfastq_template = "{output_prefix}_chunk{counter}.fastq"
            if not no_compression:
                self.outfastq_template += ".gz"
            self.next()

        def next(self):
            """Move to next chunk."""
            self.counter += 1
            logging.info("Starting chunk #{}".format(self.counter))
            self.outfastq = cancer_api.files.FastqFile.new(
                self.outfastq_template.format(**vars(self)))

        def previous(self):
            """Move to previous chunk."""
            self.counter -= 1
            logging.info("Reverting back to chunk #{}".format(self.counter))
            self.outfastq = cancer_api.files.FastqFile.open(
                self.outfastq_template.format(**vars(self)))

        def chunk_name(self):
            """Generate chunk name."""
            return "chunk{counter}".format(**vars(self))

    def split_fastq(fastq_filepath, output_dir, num_reads, num_buffer, no_compression=False):
        """Helper function for actually splitting FASTQ files.
        Simplifies the control flow depending on if one or two
        FASTQ files are specified.
        Returns intervals that were created.
        """
        # Setup
        infastq = cancer_api.files.FastqFile.open(fastq_filepath)
        current_chunk = Chunk(infastq, output_dir, no_compression)
        intervals = []
        num_accum = 0

        # The actual splitting
        for read in infastq:
            if num_accum < num_reads:
                num_accum += 1
                current_chunk.outfastq.add_obj(read)
                # Check if buffer limit reached
                if num_accum % num_buffer == 0:
                    logging.info("Flushing buffer to disk for chunk #{}".format(
                        current_chunk.counter))
                    current_chunk.outfastq.write()
            else:
                logging.info("Wrapping up chunk #{}".format(current_chunk.counter))
                current_chunk.outfastq.write()
                intervals.append(str(current_chunk.chunk_name()))
                current_chunk.next()
                num_accum = 0

        # Handle remaining reads
        # If the current_chunk is still at chunk1, just symlink
        if current_chunk.counter == 1:
            os.symlink(infastq.filepath, current_chunk.outfastq.filepath)
            intervals.append(str(current_chunk.chunk_name()))
        # Check if the number of remaining reads is above the MIN_SPILLOVER
        elif len(current_chunk.outfastq.storelist) >= MIN_SPILLOVER * num_reads:
            # If so, write out reads to current chunk
            logging.info("Writing out chunk #{}".format(current_chunk.counter))
            current_chunk.outfastq.write()
            intervals.append(str(current_chunk.chunk_name()))
        else:
            # Otherwise, write reads to previous chunk
            storelist = current_chunk.outfastq.storelist
            current_chunk.previous()
            for read in storelist:
                current_chunk.outfastq.add_obj(read)
            logging.info("Writing out chunk #{}".format(current_chunk.counter))
            current_chunk.outfastq.write()

        # Return intervals
        return intervals

    # Run split_fastq on FASTQ file(s)
    logging.info("Splitting first FASTQ file")
    intervals = split_fastq(fastq1_filepath, output_dir, args.num_reads, args.num_buffer,
                            no_compression=args.no_compression)
    if fastq2_filepath:
        logging.info("Splitting second FASTQ file")
        split_fastq(fastq2_filepath, output_dir, args.num_reads, args.num_buffer,
                    no_compression=args.no_compression)

    # Create interval_file, if applicable
    if args.interval_file:
        logging.info("Outputting intervals file")
        with open(args.interval_file, "w") as outfile:
            outfile.write("\n".join(intervals) + "\n")


if __name__ == '__main__':
    main()