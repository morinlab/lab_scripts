#~~~~~~~~~~~~
#This script takes a genelist bed and a pattern (or series of input beds containing copy number data), and generates the counts for the events from input beds that lie in the gene region. It also generates a consolidated list of these matches.
#Genelist bed is of format: 	"chr	start	end	genename"
#Input beds are of format: 	"chr	start	end	copynum"
#You can also pass in seg files (_segs.txt from Titan) but you must specify an output directory for the bed files that are generated from these as an intermediate step

#Author: Jasleen Grewal
#Date Created: March 11, 2015
#Date Modified: March 11, 2015
#~~~~~~~~~~~~

import argparse
import pandas
import re
import commands
import os

def main():
	"""Run genes_in_cn.py 
	"""

	#Specify command line arguments
	parser = argparse.ArgumentParser(description='Obtain intersect between genes and copy number events.')
	subparsers = parser.add_subparsers(dest='subcommand')

	#Subparser when input is a list of beds
	parser_files = subparsers.add_parser('in_bed')
	parser_files.add_argument('-g', '--genebed', nargs=1, type=str,required=True,help='Specify the gene bedfile')
	parser_files.add_argument('-i', '--bedfiles', nargs='+', type=str, required=True,help='Specify the pattern for the input bedfiles, or a space separated list of bedfiles')
	parser_files.add_argument('-outfile1', '--outfile1', nargs=1, type=str, required=False, default=["cn_fullcomparison_gene.txt"], help='Specify the output file for the consolidated list of matches to the genes.')
	parser_files.add_argument('-outfile2', '--outfile2', nargs=1, type=str, required=False, default=["cn_counts_gene.txt"], help='Specify the output file for the counts of matches to the genes.')

	#Subparser when input is a list of segs
	parser_files = subparsers.add_parser('in_seg')
	parser_files.add_argument('-g', '--genebed', nargs=1, type=str,required=True,help='Specify the gene bedfile')
	parser_files.add_argument('-i', '--segfiles', nargs='+', type=str, required=True, help='A space separated list of *_segs.txt files')
	parser_files.add_argument('-out', '--outdir', nargs=1, type=str, required=True, help='Output directory for the processed bed files')
	parser_files.add_argument('-outfile1', '--outfile1', nargs=1, type=str, required=False, default=["cn_fullcomparison_gene.txt"], help='Specify the output file for the consolidated list of matches to the genes.')
        parser_files.add_argument('-outfile2', '--outfile2', nargs=1, type=str, required=False, default=["cn_counts_gene.txt"], help='Specify the output file for the counts of matches to the genes.')

	#Parse command line arguments
	args = parser.parse_args()
	if args.subcommand == 'in_bed':
		input_genelist=args.genebed[0]
		num_beds = len(args.bedfiles)
		mybedlist = []
		for i in range(0,num_beds):
			mybedlist.append(args.bedfiles[i])

		fullcomparison, countscomparison = getbedintersect(input_genelist," ".join(mybedlist))
		file=open(args.outfile1[0],"w"); file.write(fullcomparison); file.close()
		file=open(args.outfile2[0],"w"); file.write(fullcomparison); file.close()
		
	if args.subcommand == 'in_seg':
		input_genelist=args.genebed[0]
		num_segs = len(args.segfiles)
		outdir = args.outdir[0]
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		outdir=re.sub("\/$","",outdir);
		mybedlist=[]
		for i in range(0,num_segs):
			segfile = args.segfiles[i]
			intrim=re.search(r"[^/]*$",segfile).group()
			intrim=re.sub("\\.txt$","",intrim); intrim=re.sub("\segs$","",intrim); intrim=re.sub("\_$","",intrim);
			outfile = "".join([outdir,"/",intrim,"_circos_bed.txt"])
			mycmd="".join(["cut -f 2-4,10 ",segfile," | sed -e 's/^/hs/' | tail -n +2 > ",outfile])
			os.system(mycmd)
			mybedlist.append(outfile)

		fullcomparison, countscomparison = getbedintersect(input_genelist," ".join(mybedlist))
                file=open(args.outfile1[0],"w"); file.write(fullcomparison); file.close()
                file=open(args.outfile2[0],"w"); file.write(fullcomparison); file.close()
		print "Your bed files are in: ", outdir

	print "The consolidated results are at: ", args.outfile1[0]
	print "The compact result counts are at: ", args.outfile2[0]
	
def getbedintersect(genebed,listofbeds):
	bedstring= listofbeds #" ".join(listofbeds)
	mycmd_consolidated = "".join(["bedtools intersect -wo -filenames -a ",genebed," -b ", listofbeds])
	mycmd_counts = "".join(["bedtools intersect -c -a ", genebed, " -b ", listofbeds])
	consolidated = commands.getoutput(mycmd_consolidated)
	counts = commands.getoutput(mycmd_counts)
#	consolidated = pybedtools.BedTool(genebed).intersect(listofbeds,wo=True,filenames=True)
#	counts = pybedtools.BedTool(genebed).intersect(listofbeds,c=True)
	return (consolidated,counts)


if __name__ == '__main__':
	main()

