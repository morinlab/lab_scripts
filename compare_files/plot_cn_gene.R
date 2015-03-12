#Pass in the consolidated results output from genes_in_cn.py

suppressWarnings(suppressMessages(library(argparser,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))
suppressWarnings(suppressMessages(library(reshape2,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))
suppressWarnings(suppressMessages(library(gplots,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))

#parser <- ArgumentParser(description='Plot the copy number states of a list of genes, for different samples. Needs output from genes_in_cn.py')
##parser$add_argument('-i',"--input", type=FileType('r'),nargs=1,help="The consolidated results file output from genes_in_cn.py")
#parser$add_argument('-out',"--output", type=FileType('w'),nargs=1,help="Destination where you wish to output the heatmap")
#args<- parser$parse_args()
#infile<- get(args$input)
#outfile<- get(args$output)

parser<- arg.parser("Plot the copy number states of a list of genes, for different samples. Needs output from genes_in_cn.py")
parser<- add.argument(parser,"--input",help="The consolidated results file output from genes_in_cn.py")
parser<- add.argument(parser,"--output",help="Destination where you wish to output the heatmap")
args<- parse.args(parser)
infile<- args$input
outfile<- args$output
print(infile)
t=read.csv(infile,sep="\t",header=TRUE)
print(dim(t))
t$cn_file=gsub(".*/.*/","",t$cn_file,perl=TRUE); t$cn_file=gsub(".txt","",t$cn_file,perl=TRUE) #Trim sample names
t$cn_state=t$cn_state-2 #Normalise CN distribution around 0 instead of around 2
tnew<-t[t$cn_file!="PT004",]; tnew<-tnew[tnew$cn_file!="PT26",] #Exclude known messy genomes

df=acast(tnew,cn_file~gene_name,value.var="cn_state") #Build a matrix from the table 

#Plot to outfile
#pdf(file=outfile)]
png(file=outfile, width=700,height=600)
heatmap.2(df,breaks=c(-3,-2,-1,0,1,2,3),key=TRUE,keysize=1,key.xlab="copy number", key.ylab="color",key.par=list(mar=c(3,4,3,3)),
          key.title="Copy Number Key (2=het)",
          key.ytickfun =function(){return(list(labels=FALSE,tick=FALSE))},
          key.xtickfun=function(){breaks <- parent.frame()$breaks
                                  return(list(at=parent.frame()$scale01(seq(from=breaks[1],to=breaks[length(breaks)],by=1))+0.1,
                                              labels=seq(from=0,to=6,by=1)))},
          col=c("dodgerblue3","dodgerblue1","white","coral1","firebrick2","firebrick4"),margins=c(10,14),
          dendrogram="col",trace="row",tracecol=rgb(1,1,1,alpha=0),linecol="black",symbreaks=FALSE)
dev.off()
