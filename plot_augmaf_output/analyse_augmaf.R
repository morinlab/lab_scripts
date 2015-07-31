# Initially written March 2015
# Last facelifted (significantly) on July 31, 2015
# By Jasleen Grewal (grewalj23@gmail.com)
# This script takes in an augmaf file, target columns, output directory 
# And plots the vaf distributions between the two samples (currently only for dual comparisons)
# May expand later if felt worth the extra effort :P
## Unfortunately currently require input in format: 
##  ("Hugo_Symbol","Chromosome","Start_Position","End_Position","Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode","t_ref_count","t_alt_count","Variant_Classification","Variant_Type")

#Req libs
suppressWarnings(suppressMessages(library(argparser,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))
suppressWarnings(suppressMessages(library(sqldf,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))
suppressWarnings(suppressMessages(library(ggplot2,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))
suppressWarnings(suppressMessages(library(grid,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))
source("multiplot.R") #Import function to plot multiple plots on same view

#Function Definitions for file parsing and processing
read_augmaf<-function(filename,sampletype,incols=selected_cols,outcols=renamed_cols){
  mydf<-read.delim(filename,stringsAsFactors=FALSE,header=TRUE)
  mydf<-mydf[incols]
  colnames(mydf) <- outcols
  mydf$vaf=mydf$t_alt_count/(mydf$t_ref_count + mydf$t_alt_count)
  mydf$sample_type = sampletype
  return(mydf)
}

compare_mafs<-function(tum_maf,plasma_maf){
  temp_tum=tum_maf[tum_maf$variant_type=="SNP",]; temp_plasma=plasma_maf[plasma_maf$variant_type=="SNP",]
  dat_overlap=sqldf("SELECT temp_tum.gene, chr, start, end,temp_tum.variant_type,temp_tum.variant_class, temp_tum.t_ref_count,temp_tum.t_alt_count, temp_tum.vaf, temp_plasma.t_ref_count,temp_plasma.t_alt_count, temp_plasma.vaf FROM temp_tum LEFT JOIN temp_plasma USING(chr,start,end)")
  colnames(dat_overlap) <- c("gene","chr","start","end","variant_type","variant_class","tum_ref","tum_alt","tum_vaf","plasma_ref","plasma_alt","plasma_vaf")
  #dat_overlap=dat_overlap[!(is.na(dat_overlap$plasma_vaf)),]
  dat_overlap[(is.na(dat_overlap$plasma_vaf)),"plasma_vaf"]=-1
  dat_overlap=dat_overlap[with(dat_overlap,order(plasma_vaf,tum_vaf)),]
  dat_overlap=dat_overlap[complete.cases(dat_overlap),]
  return(dat_overlap)
}

#Set up my argparser
parser<- arg.parser("Plot vaf distributions of paired samples from augmaf input.")
parser<- add.argument(parser,"--in_dir",help="Specify directory with augmaf files", default="/Volumes/morinlab/projects/dlbcl_normal_tumour_plasma_co/augment_maf_QCRC_out/comparative_analysis/")
parser<- add.argument(parser,"--data_dir",help="Specify directory to store comparison data files", default="/Volumes/morinlab/projects/dlbcl_relapse_exome/QCROC_cnv_snv_analyses/snv/data/")
parser<- add.argument(parser,"--plot_dir",help="Specify directory for plot files", default="/Volumes/morinlab/projects/dlbcl_relapse_exome/QCROC_cnv_snv_analyses/snv/plots/")
parser<- add.argument(parser,"--tum_min",help="Specify minimum vaf threshold. Vafs in tumour >= this value will be selected", default=0)
parser<- add.argument(parser,"--plasma_min",help="Specify minimum vaf threshold. Vafs in plasma >= this value will be selected", default=0)
parser<- add.argument(parser,"--refcov",help="Minimum number of reads in ref allele required to qualify a plasma variant (when plasma vaf=0)", default=10)
args<- parse.args(parser)

augdir= normalizePath(args$in_dir)
outputdir= normalizePath(args$data_dir)
plotdir= normalizePath(args$plot_dir)
tumour_threshold = args$tum_min
plasma_threshold = args$plasma_min
refcov = args$refcov

#Shared vars
selected_cols=c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode","t_ref_count","t_alt_count","Variant_Classification","Variant_Type")
renamed_cols=c("gene","chr","start","end","sample_type","sample_matched","t_ref_count","t_alt_count","variant_class","variant_type")

#Sample specific analysis
mysamples=c("PT003","PT011","PT012","PT018","PT32","PT39")
output_header=paste("gene","chr","start","end","variant_type","variant_class","tum_ref","tum_alt","tum_vaf","plasma_ref","plasma_alt","plasma_vaf","sample",sep="\t")
outfile_thresh=paste(outputdir,"tvaf_",tumour_threshold,"_pvaf_",plasma_threshold,".txt",sep="")
write(output_header,file=outfile_thresh,append=FALSE,sep="\t")
outfile_0=paste(outputdir,"vaf_",tumour_threshold,"_pvaf0","_cov",refcov,".txt",sep="")
write(output_header,file=outfile_0,append=FALSE,sep="\t")
outfile_na=paste(outputdir,"vaf_",tumour_threshold,"_pvafNA.txt",sep="")
write(output_header,file=outfile_0,append=FALSE,sep="\t")
#Write data
for(i in mysamples){
  sample=i
  holder_outfile=paste(outputdir,sample,"_tvaf_",tumour_threshold,"_pvaf_",plasma_threshold,".txt",sep="")
  #Read in maf files and clean them
  tumour_maffile=paste(augdir,sample,"_tumour.maf",sep="")
  ct_maffile=paste(augdir,sample,"_plasma.maf",sep="")
  tumour_maf = read_augmaf(tumour_maffile,"tumour"); 
  ct_maf=read_augmaf(ct_maffile,"plasma"); 
  #Get variants present in both 
  overlap_df=compare_mafs(tumour_maf,ct_maf)
  overlap_df$sample=sample
  write.table(overlap_df[overlap_df$tum_vaf>=tumour_threshold,],file=holder_outfile,append=FALSE,sep="\t",col.names=T, row.name=F)
  write.table(overlap_df[((overlap_df$tum_vaf>=tumour_threshold)&(overlap_df$plasma_vaf>plasma_threshold)),],file=outfile_thresh,append=TRUE,sep="\t",col.names=F,row.name=F)
  
  overlap_df_0=overlap_df[((overlap_df$tum_vaf>=tumour_threshold)&(overlap_df$plasma_ref>=refcov)&(overlap_df$plasma_vaf==0)),]
  write.table(overlap_df_0,file=outfile_0,append=TRUE,sep="\t",col.names=F, row.name=F)
  
  overlap_df_na=overlap_df[((overlap_df$tum_vaf>=tumour_threshold)&(overlap_df$plasma_vaf==-1)),]
  write.table(overlap_df_na,file=outfile_na,append=TRUE,sep="\t",col.names=F, row.name=F)
}

tumour_threshold=0.05 #vafs in tumour will be selected >= this value
plasma_threshold=0 #vafs in plasma will be selected > this value
refcov=0 #Minimum number of reads in ref allele required to qualify a plasma variant (when plasma vaf=0)
myplotarray=list()
outfile_thresh=paste(outputdir,"tvaf_",tumour_threshold,"_pvaf_",plasma_threshold,".txt",sep="")
#plot data
for(i in mysamples){
  sample=i
  #Read in maf files and clean them
  tumour_maffile=paste(augdir,sample,"_tumour.maf",sep="")
  ct_maffile=paste(augdir,sample,"_plasma.maf",sep="")
  tumour_maf = read_augmaf(tumour_maffile,"tumour"); 
  ct_maf=read_augmaf(ct_maffile,"plasma"); 
  #Get variants present in both 
  overlap_df=compare_mafs(tumour_maf,ct_maf)
  overlap_df$sample=sample
  overlap_df_0=overlap_df[((overlap_df$tum_vaf>=tumour_threshold)&(overlap_df$plasma_ref>=refcov)&(overlap_df$plasma_vaf==0)),]
  overlap_df_na=overlap_df[((overlap_df$tum_vaf>=tumour_threshold)&(overlap_df$plasma_vaf==-1)),]
  
  plot_df=overlap_df[((overlap_df$tum_vaf>=tumour_threshold)&(overlap_df$plasma_vaf>plasma_threshold)&((overlap_df$plasma_ref)>=refcov)),] #overlap_df[overlap_df$tum_vaf>=threshold,]
  #Plots!
  myplotarray[[paste(sample,"_tumour",sep="")]]=ggplot(plot_df,aes(x=tum_vaf)) + 
    geom_histogram(bindwidth=0.1,colour="cadetblue4",fill="cadetblue3") + 
    geom_density(alpha=0.3,fill="azure",colour="aquamarine") +
    ggtitle(paste(sample," Tumour vaf distribution \ntumour vaf >= ",tumour_threshold,", plasma vaf >",plasma_thresholdsep="")) +
    scale_x_continuous(breaks=seq(0,1,by=0.1)) + 
    theme(legend.position="none",text=element_text(size=20))
 myplotarray[[paste(sample,"_plasma",sep="")]]=ggplot(plot_df,aes(x=plasma_vaf)) + 
    geom_histogram(bindwidth=0.1,colour="coral4",fill="coral3") + 
    geom_density(alpha=0.3,fill="azure",colour="coral1") +
    ggtitle(paste(sample," Plasma vaf distribution \ntumour vaf >= ",tumour_threshold, ", plasma vaf >",plasma_threshold,sep="")) +
    scale_x_continuous(breaks=seq(0,1,by=0.1)) + 
    theme(legend.position="none",text=element_text(size=20))
}
jpeg(filename = paste(plotdir,"tum",tumour_threshold,"_plasma",plasma_threshold,"up_pref_",refcov,"xcov.jpeg",sep=""), width=3000,height=2500,pointsize =20, quality = 100, bg = "white", res = NA)
plotmatrix=matrix(data = 1:(ceiling(length(myplotarray)/4)*4), ncol = 4, byrow = TRUE)
plotmatrix[ceiling(length(myplotarray)/4),(4-((length(myplotarray))-(floor(length(myplotarray)/4)*4))+1):4] = 0
multiplot(plotlist=myplotarray,layout=plotmatrix)
dev.off()

#Mean vaf distributions
plot_alpha=1
mean_alpha=1
myfile=read.delim(outfile_thresh,sep="\t",stringsAsFactors=FALSE,header=TRUE)
myfile=myfile[c("tum_vaf","plasma_vaf","sample")]; colnames(myfile)=c("Tumour","Plasma","sample")
melted=melt(myfile, id.vars="sample"); colnames(melted)=c("sample","Exome_Data","value")
ggplot(melted,aes(x=sample,y=value,fill=sample,stat="sample",color=Exome_Data,alpha=plot_alpha,group=Exome_Data)) + 
  stat_summary(fun.data="median_hilow",geom="errorbar", position = position_dodge(width=0.4),width=0.2,size=0.6) +
  stat_summary(fun.y=mean,fun.ymin=mean,fun.ymax=mean, aes(colour=Exome_Data,alpha=mean_alpha,group=Exome_Data),width=2,size=1,position=position_dodge(width=0.4)) +
  ggtitle(paste(" 95% Confidence Intervals and Mean VAF Distribution in Tumour and Plasma \n For Set of Shared Variants, Tumour VAF > ",tumour_threshold,", Plasma VAF > ",plasma_threshold,sep=""))+ 
  theme(text=element_text(size=20)) +
  xlab("Sample") + ylab("Variant Allele Fraction") +
  guides(shape=FALSE,fill=FALSE,alpha=FALSE)
ggsave(file=paste(plotdir,"means_tvaf_",tumour_threshold,"_pvaf_",plasma_threshold,".jpeg",sep=""),width=14,height=9)
ggsave(file=paste(plotdir,"means_tvaf_",tumour_threshold,"_pvaf_",plasma_threshold,".pdf",sep=""),width=14,height=9)

purityfile="/Volumes/morinlab/projects/dlbcl_relapse_exome/QCROC_bed/copynum/titan_purity_comparison.txt"
purity_df=read.delim(purityfile,header=TRUE,stringsAsFactors=FALSE)
purity_df=purity_df[c("Patient","Sample","BestClust_purity","VAF_purity")]
colnames(purity_df)=c("Sample","Exome_Data","TITAN","VAF x 2")
mymelted=melt(purity_df,id.vars=c("Sample","Exome_Data"))
colnames(mymelted)=c("Sample","Exome_Data","Purity_from","Purity")
ggplot(mymelted,aes(x=Sample,y=Purity,stat="Purity_from",color=Purity_from,group=Purity_from),size=2)+
  geom_point(size=4)+facet_wrap(~Exome_Data)+ggtitle("Comparison of purity estimates from TITAN and VAFs of shared variants")+guides(alpha=FALSE)+ theme_bw()+ theme(text = element_text(size=22))
ggsave(file=paste(plotdir,"Purity_Prediction_Comparison.jpeg",sep=""),width=10,height=6)

mymelted=melt(purity_df[c("Sample","Exome_Data","TITAN")],id.vars=c("Sample","Exome_Data"))
colnames(mymelted)=c("Sample","Exome_Data","Purity_from","Purity")
ggplot(mymelted,aes(x=Sample,y=Purity,stat="Purity_from",color=Purity_from,group=Purity_from),size=2)+
  geom_point(size=4)+facet_wrap(~Exome_Data)+ggtitle("Purity estimates from TITAN")+guides(color=FALSE)
ggsave(file=paste(plotdir,"Purity_Prediction_Titan.jpeg",sep=""),width=10,height=6)



#Comparing mean vafs after and before filter
for(i in mysamples){sample=i;   tumour_maffile=paste(augdir,sample,"_tumour.maf",sep=""); 
                      ct_maffile=paste(augdir,sample,"_plasma_withplasma.maf",sep=""); 
                      tumour_maf = read_augmaf(tumour_maffile,"tumour"); 
                    ct_maf=read.delim(ct_maffile,stringsAsFactors=FALSE,header=TRUE,skip=1)
                    ct_maf<-ct_maf[selected_cols]
                    colnames(mydf) <- renamed_cols
                    ct_maf$vaf=ct_maf$t_alt_count/(ct_maf$t_ref_count + ct_maf$t_alt_count)
                    ct_maf$sample_type = "plasma"
                    
                      cat('FROM ALL: ',i, ' tumour mean vaf : ',mean(tumour_maf[tumour_maf$vaf>0,"vaf"],na.rm=TRUE), '; plasma mean vaf : ',mean(ct_maf[ct_maf$vaf>0,"vaf"],na.rm=TRUE),"\n"); 
                      cat(i, ' tumour median vaf : ',mean(myfile[myfile$sample==sample,"tum_vaf"]), '; plasma median vaf : ',mean(myfile[myfile$sample==sample,"plasma_vaf"]),"\n"); }
#----plot 1
names=c("Tumour","Plasma")
data=data.frame()
ggplot()
myplot=ggplot(mydf,aes(x=cn,y=reads_per_base,fill=sample,stat="sample",color=sample,alpha=plot_alpha)) + stat_summary(fun.data="median_hilow",geom="pointrange", position = position_dodge(width=0.4)) 
myplot + facet_wrap(~sample_type, scales="free") +  stat_summary(fun.y=mean, aes(colour=sample,alpha=mean_alpha,group=sample),size=1,geom="line")
ggsave(file="plots/nonnormalized_pertype_rd.pdf",width=14,height=9)

