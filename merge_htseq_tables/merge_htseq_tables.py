#!/usr/bin/env python
#####################################################
#                                                   #
# Author: Dave Wheeler                              #
#                                                   #
# Purpose: merge count tables                       #
#                                                   #
# Usage: python merge_tables guide_file             #
#####################################################
#looks for text guide_file
#that contains files to be merged with column headers (space separted) 
#ie
#file1.counts untreated1
#file2.counts untreated2
#file3.counts treated1
#file4.counts treated2
#this will generated a tab separated table like this
#
#gene untreated1 untreated2 treated1 treated2
#gene1 0 0 0 0 
#gene2 1 0 11 10
#.......
##############################################
import sys
try:
	infile = open(sys.argv[1])
except IndexError:
	print "No guide file provided"
	sys.exit()
	
#make dict of genes with list of counts
#list is ordered so treatments will be preserved.
#genes = {'gene1':[1,2,3,4]}
#header keeps track of treatment order, will be as read from config
col_header = []
genes = {}	

outfile = open('merged_counts.txt','w')

for line in infile:
	filename,header = line.strip().split(' ')
	try:
		data_f = open(filename)
	except IOError:
		print "%s can't be found?"%filename
		sys.exit()
		
	col_header.append(header)	
	
	#read file and add gene and counts to the dict
	for line in data_f:
		gene,count = line.strip().split('\t')
		if gene not in genes:
			genes[gene] = [count]
		else:
			genes[gene].append(count)
	#important to close file
	data_f.close()	
	
infile.close()

outfile.write('gene\t'+'\t'.join(col_header)+'\n')

for gene in genes:
	data = genes[gene]
	#make sure each treatment has a count for this gene
	#this should catch most errors
	try:
		assert len(data) == len(col_header)
	except AssertionError:
		print "one of the treatment or genes is missing or extra"
		print "data, found the problem here:"
		print gene,data
		print "while %s columns of treatments given" %len(col_header)
		sys.exit()
		
	out_data = gene+'\t'+'\t'.join(data)+'\n'
	outfile.write(out_data)
	
outfile.close()
print "Merged table is 'merged_counts.txt'"	
