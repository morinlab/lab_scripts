#example usage: python genes_in_cn.py in_bed -g circos-0.67-4/data/circos_label/genes/list_genes.txt -i circos-0.67-4/data_pan_qcroc/circos_data/*.txt -oprefix allpts
#Function defs
sample_format<-function(mydf_col){ #File output when you go readintersect(cvexons,sample[chr,start,end,reads])
  mydf_col=gsub(".*/.*/","",mydf_col,perl=TRUE); mydf_col=gsub(".txt","",mydf_col,perl=TRUE) #Trim sample names
  mydf_col=gsub("PAN","CH-",mydf_col,perl=TRUE); mydf_col=gsub("PT","QC2-",mydf_col,perl=TRUE) #Trim sample names
  mydf_col=gsub("QC2-00(?=.$)","QC2-0",mydf_col,perl=TRUE); mydf_col=gsub("QC2-0(?=..$)","QC2-",mydf_col,perl=TRUE);
  mydf_col=gsub("CH-(?=..$)","CH-0",mydf_col,perl=TRUE);mydf_col=gsub("^(?=...$)","MT-",mydf_col,perl=TRUE);
  return(mydf_col);
}

matchdf<-function(sortMe,sortBy){
  sorted=sortBy
  for(i in 1:ncol(sortBy)){
    gene=colnames(sortBy)[i]
    for(j in 1:nrow(sortBy)){
      sample=rownames(sortBy)[j]
      if((sample %in% rownames(sortMe))&(gene %in% colnames(sortMe))){
        sorted[j,i]=sortMe[sample,gene]
      }else{
        sorted[j,i]=""
      }}}
  return(sorted)
}
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
cat("MAKING R PLOT \nReading from.....", infile,"\n")
t=read.csv(infile,sep="\t",header=TRUE)
t$cn_file=gsub(".*/.*/","",t$cn_file,perl=TRUE); t$cn_file=gsub(".txt","",t$cn_file,perl=TRUE) #Trim sample names
t<-t[t$cn_file!="PT004",]; t<-t[t$cn_file!="PT26",] #Exclude known messy genomes
t$cn_file=sample_format((t$cn_file))
t$cn_state=t$cn_state-2 #Normalise CN distribution around 0 instead of around 2
tnew<-t
df=acast(tnew,cn_file~gene_name,value.var="cn_state",fun.aggregate=mean,fill=0) #Build a matrix from the table 

df2 <- (as.matrix(read.table("reference_sv_data/snv_indel_data.txt", header=TRUE, sep = "\t",row.names = 1,as.is=TRUE)))
df2<-ifelse(df2==1,"\\",ifelse(df2==2,"x",ifelse(df2==3,"*",""))) #1 is snv,  2 is indel, 3 is snv and indel

rownames(df2) =sample_format(rownames(df2))
df2.ordered<-matchdf(df2,df)

#Plot to outfile
#pdf(file=outfile)]
png(file=paste(outfile,".png",sep=""), width=1500,height=1000,res=100)
heatmap.2(df,breaks=c(-3,-2,-1,0,1,2,3),key=TRUE,keysize=1,key.xlab="copy number", key.ylab="color",key.par=list(mar=c(3,4,3,3)),
          key.title="Copy Number Key (2=het)",
          key.ytickfun =function(){return(list(labels=FALSE,tick=FALSE))},
          key.xtickfun=function(){breaks <- parent.frame()$breaks
                                  return(list(at=parent.frame()$scale01(seq(from=breaks[1],to=breaks[length(breaks)],by=1))+0.1,
                                              labels=seq(from=0,to=6,by=1)))},
          col=c("dodgerblue3","dodgerblue1","white","coral1","firebrick2","firebrick4"),margins=c(10,14),
          dendrogram="col",trace="none",tracecol=rgb(1,1,1,alpha=0),colsep=1:ncol(df),rowsep=1:nrow(df),sepcolor="grey",sepwidth=c(0.01,0.01),symbreaks=FALSE,par(cex.main=1.5),
          linecol="black")
          #trace="both",tracecol=rgb(1,1,1,alpha=0),linecol="black",symbreaks=FALSE,par(cex.main=1.5))
dev.off()

png(file=paste(outfile,"_fullclust.png",sep=""), width=1500,height=1000,res=100)
heatmap.2(df,breaks=c(-3,-2,-1,0,1,2,3),key=TRUE,keysize=1,key.xlab="copy number", key.ylab="color",#key.par=list(mar=c(3,4,3,3)),
          key.title="Copy Number Key (2=het)",
          key.ytickfun =function(){return(list(labels=FALSE,tick=FALSE))},
          key.xtickfun=function(){breaks <- parent.frame()$breaks
                                  return(list(at=parent.frame()$scale01(seq(from=breaks[1],to=breaks[length(breaks)],by=1))+0.1,
                                              labels=seq(from=0,to=6,by=1)))},
          col=c("dodgerblue3","dodgerblue1","white","coral1","firebrick2","firebrick4"),
          dendrogram="both",trace="none",tracecol=rgb(1,1,1,alpha=0),
          colsep=1:ncol(df),rowsep=1:nrow(df),sepcolor="grey",sepwidth=c(0.01,0.01),symbreaks=FALSE)
dev.off()

png(file=paste(outfile,"_sv.png",sep=""), width=1500,height=1000,res=100)
heatmap.2(df,cellnote=df2.ordered,notecol="black",notecex=1.3,na.color=par("bg"),breaks=c(-3,-2,-1,0,1,2,3),key=TRUE,keysize=1,key.xlab="copy number", key.ylab="color",key.par=list(mar=c(3,4,3,3)),
          key.title="Copy Number Key (2=het)",
          key.ytickfun =function(){return(list(labels=FALSE,tick=FALSE))},
          key.xtickfun=function(){breaks <- parent.frame()$breaks
                                  return(list(at=parent.frame()$scale01(seq(from=breaks[1],to=breaks[length(breaks)],by=1))+0.1,
                                              labels=seq(from=0,to=6,by=1)))},
          col=c("dodgerblue3","dodgerblue1","white","coral1","firebrick2","firebrick4"),margins=c(10,14),
          dendrogram="both",trace="none",tracecol=rgb(1,1,1,alpha=0),
          colsep=1:ncol(df),rowsep=1:nrow(df),sepcolor="grey",sepwidth=c(0.01,0.01),symbreaks=FALSE)
dev.off()

cat("Note legend for ",paste(outfile,"_sv.png",sep="")," is as follows: \n \t \\ = snv, x = indel, * = snv + indel \n")

