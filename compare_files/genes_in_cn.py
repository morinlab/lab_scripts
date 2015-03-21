#~~~~~~~~~~~~
#This script takes a genelist bed and a pattern (or series of input beds containing copy number data), and generates the counts for the events from input beds that lie in the gene region. It also generates a consolidated list of these matches.
#Genelist bed is of format: 	"chr	start	end	genename"
#Input beds are of format: 	"chr	start	end	copynum"
#You can also pass in seg files (_segs.txt from Titan) but you must specify an output directory for the bed files that are generated from these as an intermediate step
#Note that using in_seg will only generate your bed files. You must subsequently pass these bed files to this script with in_bed, in order to plot the copy number distribution.
#Author: Jasleen Grewal
#Date Created: March 11, 2015
#Date Modified: March 20, 2015
#~~~~~~~~~~~~

import argparse
import pandas
import re
import commands
import subprocess
import os

plotscript = "plot_cn_gene.R"
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
	parser_files.add_argument('-oprefix', '--outfile_prefix', nargs=1, type=str, required=False, default=["cn_gene"], help='Specify the prefix for the output files.')

	#Subparser when input is a list of segs
	parser_files = subparsers.add_parser('in_seg')
    #parser_files.add_argument('-g', '--genebed', nargs=1, type=str,required=True,help='Specify the gene bedfile')
	parser_files.add_argument('-i', '--segfiles', nargs='+', type=str, required=True, help='A space separated list of *_segs.txt files')
	parser_files.add_argument('-out', '--outdir', nargs=1, type=str, required=True, help='Output directory for the processed bed files')
        #parser_files.add_argument('-oprefix', '--outfile_prefix', nargs=1, type=str, required=False, default=["cn_gene"], help='Specify the prefix for the output files.')

	#Parse command line arguments
	args = parser.parse_args()
	input_bedlist=""
	filecount=0

	if args.subcommand == 'in_bed':
		input_genelist=args.genebed[0]
		output_prefix=args.outfile_prefix[0]
		out_consolidated="".join([output_prefix,"_consolidated.txt"])
		out_counts="".join([output_prefix,"_counts.txt"])
		out_image=output_prefix

		filecount = len(args.bedfiles)
		mybedlist = []
		for i in range(0,filecount):
			mybedlist.append(args.bedfiles[i])
		input_bedlist= " ".join(mybedlist)
		fullcomparsion, countscomparison = getbedintersect(input_genelist, input_bedlist,out_consolidated,out_counts,out_image,filecount)
		print "The consolidated results are at: ", out_consolidated
		print "The compact result counts are at: ", out_counts
		print "View it graphically! : ", "".join([out_image,"_*.png"])

		
	if args.subcommand == 'in_seg':
		filecount = len(args.segfiles)
		outdir = args.outdir[0]
		if not os.path.exists(outdir):
	            os.makedirs(outdir)
		outdir=re.sub("\/$","",outdir);
		mybedlist=[]
		for i in range(0,filecount):
			segfile = args.segfiles[i]
			intrim=re.search(r"[^/]*$",segfile).group()
			intrim=re.sub("\\.txt$","",intrim); intrim=re.sub("\segs$","",intrim); intrim=re.sub("\_$","",intrim);
			outfile = "".join([outdir,"/",intrim,"_circos_bed.txt"])
			mycmd="".join(["cut -f 2-4,10 ",segfile," | sed -e 's/^/hs/' | tail -n +2 > ",outfile])
			os.system(mycmd)
			mybedlist.append(outfile)
		input_bedlist= " ".join(mybedlist)
 		print "Your bed files are in: ", outdir


	
def getbedintersect(genebed,listofbeds,out1,out2,outimg,numfiles):
	bedstring= listofbeds #" ".join(listofbeds)
	mycmd_consolidated = "".join(["bedtools intersect -wo -filenames -a ",genebed," -b ", listofbeds])
	mycmd_counts = "".join(["bedtools intersect -c -a ", genebed, " -b ", listofbeds])
	mycmd_image = "".join(["Rscript ",plotscript," -i ",out1," -o ", outimg])
	#Run commands to generate data files
	consolidated = commands.getoutput(mycmd_consolidated)
	counts = commands.getoutput(mycmd_counts)
	#Write to file
	file=open(out1,"w"); file.write("\t".join(["gene_chr","gene_start","gene_end","gene_name","cn_file","cn_chr","cn_start","cn_end","cn_state","bases_overlap"])); file.write("\n"); file.write(consolidated); file.close()
	file=open(out2,"w"); file.write("\t".join(["gene_chr","gene_start","gene_end","gene_name","number_of_matched_samples"])); file.write("\n"); file.write(counts); file.write("".join(["\nTotal files are ",str(numfiles),"\n"])); file.close()
	#Generate plot
	subprocess.call(mycmd_image,shell=True)
	return (consolidated,counts)


if __name__ == '__main__':
	main()

