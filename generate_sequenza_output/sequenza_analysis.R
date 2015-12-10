
#Date modified: 25th August, 2015
#Author: Jasleen Grewal (jkg21@sfu.ca)
#Import Libraries
if(!require("sequenza"))install.packages("sequenza",repos="http://cran.rstudio.com/")
if(!require("argparse"))install.packages("argparse",repos="http://cran.rstudio.com/")
if(!require("copynumber"))(source("http://bioconductor.org/biocLite.R"))
if(!require("copynumber"))(biocLite("copynumber"))

suppressWarnings(suppressMessages(library(argparse,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))
suppressWarnings(suppressMessages(library(sequenza,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))


#Specify command line arguments
parser <- ArgumentParser()
parser$add_argument("-c", "--chr", action="store_false", dest="chr_list",help="Analyze per chromosome.")
parser$add_argument("--input",help="Specify the .seqz.gz input.")
parser$add_argument("--output",help="Destination directory for analysis results.")
parser$add_argument("--sample_id",help="Specify Sample ID for this input.")
parser$add_argument("--gender", help="Specify patient's gender. f or m")
parser$add_argument("--min_reads",default=40,type="integer",help="Specify min reads in tumour bam")
parser$add_argument("--min_normal_reads",default=10,type="integer",help="Specify min reads in normal bam")
parser$add_argument("--cellularity_limit",default=1,type="double",help="Specify x in iterations of cellularity from 0.1 to x, step size 0.01. [Default \"%(default)s\"]")
parser$add_argument("--ploidy_limit",default=7,type="integer",help="Specify x in iterations of ploidy from 1 to x, step size 0.1. [Default \"%(default)s\"]")
parser$add_argument("--max_cn",default=20,type="integer",help="Specify max copy number to consider in the model. [Default \"%(default)s\"]")

#Process arguments
args <- parser$parse_args()
infile<- args$input
outdir<- args$output
sample_id<- args$sample_id

if(args$chr_list){
  chromosome_list=seq(1,22,1)
  chromosome_list=unlist(list(chromosome_list,"X","Y"))
}else{
  chromosome_list=NULL
}

#Analysis step 1 - read raw data, normalize depth ratio for GC content bias, 
#perform allele-specific segmentation, filter noisy mutations, bin raw data for plotting.
#RETURNS a single list object.
my_seq_extract = sequenza.extract(
    file = infile
)

#Analysis step 2 - calculate log posterior probability for all pairs of candidate ploidy
#and cellularity parameters
my_seq_fit = sequenza.fit(
    my_seq_extract
)
#Prepare inputs for Sequenza (Python Pre Processing Script)

#Analysis step 3 - save output objects in output directory
sequenza.results(
    my_seq_extract, 
    my_seq_fit, 
    out.dir = outdir, 
    sample.id = sample_id
)
