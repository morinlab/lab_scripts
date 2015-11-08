#!/usr/bin/env python

"""
merge_mpileup.py
==========
Description
Merges mpileup (.pileup) files in a given directory.

Dependencies
-None

Inputs
-Directory containing pileup files.

Outputs
-Single merged pileup file.

Known Issues
-None

"""

import os.path
import glob
import subprocess
import argparse
import re


def main():
	args = parse_args()
	input_dir = args.input_dir
	output_file = args.output_file
	
	if not os.path.exists(input_dir):
		raise ValueError("{} does not exist.".format(input_dir))
	
	if not os.path.isdir(input_dir):
		raise ValueError("{} is not a directory.".format(input_dir))
	
	pileups = glob.glob(os.path.join(input_dir,"*.pileup"))
        if not len(pileups):
                raise ValueError("No pileup files found in {}".format(input_dir))

	order = ['1', '2', '3', '4', '5', '6', '7', 'X', '8', '9', '10', 
                 '11', '12', '13', '14', '15', '16', '17', '18', '20', 
                 'Y', '19', '22', '21']
	ordered_pileups = []
	for i in order:
		s = ".+_{}.pileup".format(i)
		p = re.compile(s)
		ordered_pileups.extend([m.group(0) for l in pileups for m in [p.search(l)] if m])
	
	with open(output_file, 'w') as outfile:
		for pileup in ordered_pileups:
			with open(pileup) as infile:
				for line in infile:
					outfile.write(line)
	return

def parse_args():
	parser = argparse.ArgumentParser()
	
	parser.add_argument("input_dir", help="Directory containing mpileup files.")
	parser.add_argument("output_file", help="Merged mpileup file.")
	# Add option to specify file extension for pileup files, default to
	# '.pileup'
	return parser.parse_args()

	
if __name__ == '__main__':
	main()
