#~~~~~~~~~
#This script takes in two files and outputs the intersecting entries, as well as entries unique to each input file, in three separate output files.
#Assumptions about input data: 
	#If comparing two files based on chr, start, end, then both should have columns: "Chromosome", "Start_Position", "End_Position"
	#If filtering one variant file based on a list of genes, variant file has 'Hugo_Symbol' column for filtering, and the first column in input file for list of genes has hugo symbols. 
	#Assumes no header for gene list file.
#You can specify which row to take as the header in each of the two inputs (default being 0,0). 
#Tested on these paired inputs : {(augmaf,augmaf),(augmaf,genelist)}

#Author: Jasleen Grewal.
#Date Created: February 27th, 2015.
#Date Last Modified: April 16th, 2015.
#Sample usage: python filecompare.py filecompare -i file1.maf file2.maf -out file1ID file2ID -m -header 0 1
#Sample usage: python filecompare.py genefilter -file file.maf -genes genelist.txt -out fileID -header 1

#Output from 'filecompare': file1ID_only.maf (variants unique to file1), file2ID_only.maf (variants unique to file2), comparison_shared.maf (variants common to both inputs), comparison_merged (Remove duplicate rows by criteria : Chromosome-Start Position-End Position) 
#Output from 'genefilter': fileID_filtered.maf (variants for genes from passed list), fileID_exclude_genes.maf (variants for genes not in the passed list).
#~~~~~~~~~~

import pandas
import argparse
import numpy

def main():
	"""Run augmentmaf_filter.py
	"""
	# Specify command line arguments
	parser = argparse.ArgumentParser(description='Obtain the intersect between two files.')
	subparsers = parser.add_subparsers(dest='subcommand')

	#subparser for 2 file comparison
	parser_files = subparsers.add_parser('filecompare')
	parser_files.add_argument('-i', '--input', nargs=2, type=argparse.FileType('r'),required=True,
                        help='Specify the two files whose intersect you\'re interested in.')
	parser_files.add_argument('-out', '--output_file', nargs=2, type=str, required=True,help='Specify the prefixes to identify output files for variants unique to file1 and file2 respectively.')
	parser_files.add_argument('-m', '--merge', action='store_true', default=False, help='Instead of finding the intersecting rows, this script merges both input files such that there are no duplicate rows.')
	parser_files.add_argument('-header', '--header_row',nargs=2, type=int, default=[0,0],help='If your input files have different header indices, pass them here.')


	#subparser for 1 infile and 1 genelist filter
	parser_filter = subparsers.add_parser('genefilter')
	parser_filter.add_argument('-i', '--input', nargs=1, type=argparse.FileType('r'),required=True, help="Specify the file you want to filter")
	parser_filter.add_argument('-genes', '--genelist', nargs=1, type=argparse.FileType('r'), required=False, default = ['/home/jgrewal/scripts/ref/genes_lymphoma.txt'],help="Specify the gene list you wish to filter against")
	parser_filter.add_argument('-out', '--output_file', nargs=1, type=str, required=True, help='Specify the prefix to identify output file.')
	parser_filter.add_argument('-header', '--header_row', nargs=1, type=int, default=[0], help='If your input file has a different header index, pass it here')

    	# Parse command line arguments
	args = parser.parse_args()
	if args.subcommand == 'filecompare': 
		input_maf_1 = pandas.DataFrame.from_csv(args.input[0], sep="\t",index_col=None,header=args.header_row[0]);
		input_maf_1 = updatecoltypes(input_maf_1);
		input_maf_2 = pandas.DataFrame.from_csv(args.input[1], sep="\t",index_col=None,header=args.header_row[1]);
		input_maf_2 = updatecoltypes(input_maf_2);
		is_merging = args.merge
	
		#Get your intersects and unique values
		maf1only = onlyInLeft(input_maf_1,input_maf_2)
		maf2only = onlyInLeft(input_maf_2,input_maf_1)
		intersect = InBoth(input_maf_1,maf1only)

		maf1only.to_csv((args.output_file[0]+"_only.maf"), sep="\t", index=False)
		maf2only.to_csv((args.output_file[1]+"_only.maf"),sep="\t", index=False)
		intersect.to_csv(("comparison_shared.maf"),sep="\t", index=False)

		if (is_merging):
			mergefile = pandas.concat([maf1only,maf2only]).drop_duplicates().reset_index(drop=True)
			mergefile.to_csv(("comparison_merged.maf"),sep="\t", index=False)

	if args.subcommand == 'genefilter':
		input_maf = pandas.DataFrame.from_csv(args.input[0], sep="\t",index_col=None,header=args.header_row[0])
		input_list = pandas.DataFrame.from_csv(args.genelist[0], sep="\t", index_col=None, header=None)
		filtered = GeneFilter(input_maf,input_list)
		excludegenes = GeneFilter_exclude(input_maf,input_list)
		filtered.to_csv((args.output_file[0] + "_filtered.maf"), sep="\t",index=False)
		excludegenes.to_csv((args.output_file[0] + "_exclude_genes.maf"),sep="\t",index=False)

def onlyInLeft(table1,table2):
        return table1[~((table1.Chromosome.isin(table2.Chromosome)) & (table1.Start_Position.isin(table2.Start_Position)) & (table1.End_Position.isin(table2.End_Position)))]

def InBoth(table1,table2):
# return table1[((table1.Chromosome.isin(table2.Chromosome)) & (table1.Start_Position.isin(table2.Start_Position)) & (table1.End_Position.isin(table2.End_Position)))]
    tempmerge = pandas.concat([table1,table2]).reset_index(drop=True)
    return tempmerge.groupby(["Chromosome","End_Position","Start_Position"]).filter(lambda x: len(x) ==1)

def GeneFilter(table1,table2):
	return table1[((table1.Hugo_Symbol.isin(table2.loc[0:,0])))]

def GeneFilter_exclude(table1,table2):
        return table1[~((table1.Hugo_Symbol.isin(table2.loc[0:,0])))]

def updatecoltypes(mytable):
    temptable=makecolint(mytable,["Start_Position","End_Position"])
    temptable=makecolstr(temptable,["Chromosome"])
    return temptable

def makecolint(table1,listofcols):
        table1[listofcols] = table1[listofcols].astype(numpy.int64)
        return table1

def makecolstr(table1,listofcols):
    table1[listofcols] = table1[listofcols].astype(str)
    return table1

if __name__ == '__main__':
	main()
