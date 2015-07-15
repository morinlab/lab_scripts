#This script takes in a directory of log files
#Greps filenames of format <taskname>.o<job_id>
#And queries runtime and additional stats for each job id
#Output is <taskname>.qacct<job_id> in the provided directory
#Ofcourse, you must have a queue manager in the cluster for this to work :P

#Author: Jasleen Grewal (grewalj23@gmail.com)
#Date Created: July 15, 2015
#Date Last Modified : July 15, 2015

import os
import argparse
import re
from os import stat
from pwd import getpwuid
import glob
import fnmatch
from subprocess import call

def main():
	"""Run get_job_stats.py
	"""
	# Specify command line arguments
	parser = argparse.ArgumentParser(description='Generate .qacct files for each job')
	parser.add_argument('-idir','--indir', nargs=1, action='store',help='Log directory with list of tasks')

	args = parser.parse_args()
	#Get list of files
	my_logs_dir = args.indir[0]
	
	for file in os.listdir(my_logs_dir):
		if fnmatch.fnmatch(file,'*.o*'):
			file_jobid = re.search('(?<=\.o)[0-9]*',file).group(0)
			file_name = re.search('.*(?=\.o)',file).group(0)
			file_owner = find_owner(my_logs_dir+'/'+file)
			my_out_file=open(my_logs_dir+"/"+file_name+".qacct"+file_jobid,"w")
			return_call = call("qacct -j "+file_jobid+" -o " + file_owner,shell=True,stdout=my_out_file)			

def find_owner(myfile):
	return getpwuid(stat(myfile).st_uid).pw_name

if __name__ == '__main__':
	main()
