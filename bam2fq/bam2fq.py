#!/usr/bin/env python

"""
bam2fq.py
==========
Description
Converts each read in a SAM file into FASTQ format while pairing
together mate end pair reads. Also shuffles the reads to prevent
bias during alignment.

Dependencies
-None

Inputs
-SAM-formatted file of read mappings

Outputs
-FASTQ file for paired reads
-FASTQ file for unpaired reads

Known Issues
-Refactor code to adhere to DRY
-FASTQ reads will not be shuffled

"""

import argparse
import random


def main():
	"""Parse BAM file and convert to FASTQ"""

	# Argument parsing
	args = parse_args()
	infile_bam = args.bam[0]
	unpaired_outfile = args.unpaired_out[0]
	paired_outfile = args.paired_out[0]
	# Initalize data structures
	read_pair_dict = {}
	unpaired_reads = {}
	paired_read_names = []
	# Look over each line in bam
	for line in infile_bam:
		bam_line = line.split("\t")
		qname = bam_line[0]
		sam_flag_list = decompose_flag(bam_line[1], [1, 16, 64, 128])
		seq = bam_line[9]
		quality = bam_line[10]

		# Reverse complement sequence and reverse quality string
		# if revcomp flag is set
		if sam_flag_list[1]:
			seq = revcomp(seq)
			quality = quality[::-1]

		# if read is unpaired, store in unpaired_reads dict
		# do I append a read pair number (i.e., \1 or \2)?
		if not sam_flag_list[0]:
			unpaired_reads[qname] = (seq, quality)
			continue

		# Add qname as key if it doesn't exist.
		# If it already exists, means the current read is part of
		# a pair.
		if qname not in read_pair_dict.keys():
			read_pair_dict[qname] = ([], [])
		else:
			paired_read_names.append(qname)

		# Add sequence and quality to array corresponding to the
		# read
		if sam_flag_list[2]:
			sort_reads(read_pair_dict[qname], 0, seq, quality)
		elif sam_flag_list[3]:
			sort_reads(read_pair_dict[qname], 1, seq, quality)

	else:
		infile_bam.close()

	# Identifies reads with only forward or reverse read
	for unpaired in set(read_pair_dict.keys()) - set(paired_read_names):
		if len(read_pair_dict[unpaired][0]):
			unpaired_reads[unpaired] = (x[0][0], x[0][1])
		else:
			unpaired_reads[unpaired] = (x[1][0], x[1][1])

	# Shuffle paired reads prior to output
	random.shuffle(paired_read_names)

	# Output paired reads to file
	for k in paired_read_names:
		t = [k, read_pair_dict[k][0][0], read_pair_dict[k][0][1], read_pair_dict[k][1][0], read_pair_dict[k][1][1]]
		outstring = '@{0}/1\n{1}\n+\n{2}\n@{0}/2\n{3}\n+\n{4}\n'
		paired_outfile.write(outstring.format(*t))
	else:
		paired_outfile.close()

	# Output unpaired reads to file
	for k in unpaired_reads.keys():
		t = [k, unpaired_reads[k][0], unpaired_reads[k][1]]
		unpaired_outfile.write('@{0}\n{1}\n+\n{2}\n'.format(*t))
	else:
		unpaired_outfile.close()


# Identifies the bits that are set in the SAM flag
def decompose_flag(flag, bits=[2**x for x in range(12)]):
	flag_list = [bool(int(flag) & x) for x in bits]
	return flag_list


# Reverse complements a sequence
def revcomp(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	rev = "".join(complement.get(base, base) for base in reversed(seq))
	return rev


# Appends items to a list at index 'i' in a given tuple reference
# ***No return value***
def sort_reads(tuple_to_append, list_index, seq, quality):
	tuple_to_append[list_index].append(seq)
	tuple_to_append[list_index].append(quality)


def parse_args():
	"""Parse command-line arguments"""
	# Setup
	parser = argparse.ArgumentParser()
	# Positional arguments
	parser.add_argument('bam', nargs=1, type=argparse.FileType('r'), help = 'Specify paired end BAM file.')
	parser.add_argument("paired_out", nargs=1, type=argparse.FileType("w"))
	parser.add_argument("unpaired_out", nargs=1, type=argparse.FileType("w"))
	# Parsing
	args = parser.parse_args()
	return args


if __name__ == '__main__':
	main()
