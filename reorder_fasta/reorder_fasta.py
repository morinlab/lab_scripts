#!/usr/bin/env python

import argparse
import os
import tempfile

# Written by Christopher Rushton
# Please forward all cookies/comments to ckrushto@sfu.ca


def is_valid_file(filepath, parser):

    if os.path.exists(filepath):
        return filepath
    else:
        raise parser.error("Unable to locate \'%s\': No such file or directory" % filepath)


def get_args():
    """
    Process command line arguments
    """

    parser = argparse.ArgumentParser(description="Reorders a FASTA file based on a user-specified order")
    parser.add_argument("-f", "--fasta", required=True, type=lambda x: is_valid_file(x, parser), help="Input FASTA file (to be reordered)")
    parser.add_argument("--order", required=True, type=lambda x: is_valid_file(x, parser), help="An input txt file specifying the sort order for FASTA entries")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output FASTA file")
    parser.add_argument("--tmp", type=str, default=None, help="Temp folder to use when re-ordering FASTA [Default: Create a system temp folder]")
    return parser.parse_args()


def split_fasta(fasta_file, tmp_dir):
    """
    Splits a FASTA file by contig, creating an output file for each contig in the input FASTA

    :param fasta_file: A string containing a filepath to the FASTA file
    :param tmp_dir: A string containing a filepath to a temporary directory, which must exist
    """

    first_line = True
    outfile = None
    contig_name = None
    fasta_file_names = {}
    with open(fasta_file) as f:
        for line in f:

            line = line.rstrip("\n").rstrip("\r")
            # Sanity check input file format
            if first_line:
                first_line = False
                if not line.startswith(">"):
                    raise AttributeError("Input FASTA file \'%s\' does not appear to be in FASTA format (does not start with \">\" in the first line)" % fasta_file)
            if line.startswith(">"):
                # New FASTA entry. Create a new file for this entry
                if outfile is not None:
                    # Close the existing FASTA. We have finished processing it
                    outfile.close()
                # Use the name of this entry/contig as the fasta name
                fasta_name = line.split()[0]
                fasta_name = fasta_name.replace(">", "")
                fasta_file_name = tmp_dir + os.sep + fasta_name + ".fa"
                fasta_file_names[fasta_name] = fasta_file_name
                outfile = open(fasta_file_name, "w")

            # Write FASTA contents to the appropriate file
            outfile.write(line)
            outfile.write(os.linesep)

    return fasta_file_names


def main(args=None):

    if args is None:
        args = get_args()


    # Prepare temp dir for storing split FASTAs
    if args.tmp is None:
        tmp_dir = tempfile.mkdtemp()
    else:
        # Don't overwrite/use existing directories
        if os.path.exists(args.tmp):
            raise AttributeError("Specified temp directory already exists: \'%s\'" % args.tmp)
        tmp_dir = args.tmp
        os.mkdir(args.tmp, mode=0o700)

    # Step 1: Split the FASTA file by entry
    contig_files = split_fasta(args.fasta, tmp_dir)

    # Step 2: Load contig order
    contig_order = []
    with open(args.order) as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            # Sanity check this contig exists in the reference
            if not line in contig_files:
                raise AttributeError("FASTA entry \'%s\' was specified in the --order file, but is not present in the FASTA file" % line)
            contig_order.append(line)

    # Step 3: Write out contigs in the specified order
    with open(args.output, "w") as o:
        for contig in contig_order:

            # Get the appropriate FASTA file for this entry, and write it out
            fasta_file = contig_files[contig]
            with open(fasta_file) as f:
                for line in f:
                    o.write(line)
            # For my own personal sanity, delete the FASTA contig file after reading it for the
            # output file
            os.remove(fasta_file)

    # Finally, clean up the temp folder
    os.rmdir(tmp_dir)


if __name__ == "__main__":
    main()
