#This file takes in input file list and file type of list members
#And outputs the sum of the counts depending on type of file

import commands
import os
import re
import argparse
import subprocess

def cmd_readcount(infile, filetype, outfile, whichsam):
	if filetype == "fastq":
        	mycmd = "cat " + infile + " | wc -l | xargs -n 1 bash -c 'echo $(($1/4))' args " + " >> " + outfile
	elif filetype == "bam":
		mycmd = whichsam + " view -F 2048 -c " + infile + " >> " + outfile
	elif filetype == "fq.gz":
		mycmd = "zcat -f " + infile + " | wc -l | xargs -n 1 bash -c 'echo $(($1/4))' args " + " >> " + outfile
	return mycmd

def main():
	"""Run count_comparefiles.py
	"""

	#Specify command line arguments
	parser = argparse.ArgumentParser(description="Get total counts for file list, depending on file type")
        parser.add_argument('-t', '--input_type', nargs=1, type=str, required=True, help="One of fastq or bam or fq.gz")
        parser.add_argument('-o', '--countfile', nargs=1, type=str, required=True, help="Name of output file with one count per file list item")
	parser.add_argument('-i', '--input_list', nargs='+', type=str, required=True, help="List of input file names")
	parser.add_argument('--samtools', nargs=1, type=str, required=False, default=["samtools"], help="Path to samtools exec")

	#Parse command line arguments
	args = parser.parse_args()
	count_sums = 0
	input_files = args.input_list
	my_file_type = args.input_type[0]
	myoutfile = args.countfile[0]
	runcmd = []
	mysamtools = args.samtools[0]

	for myinfile in args.input_list:
		runcmd.append(cmd_readcount(myinfile,my_file_type,myoutfile,mysamtools))
	
	myoutputcmds = " ; ".join(runcmd)
	print myoutputcmds
	os.system(myoutputcmds)
	return myoutputcmds

if __name__ == '__main__':
	main()
