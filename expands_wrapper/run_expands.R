#!/usr/bin/env Rscript
#simple wrapper to run the expands software using Sequenza-flavour copy number results as input
#there is no error checking and the script blindly assumes you have the appropriate dependencies installed in the R on your path and will run it as follows:
# this_script.R sequenza_or_titan_segments_file.txt augmented_maf_file.maf sample_name max_score precision input_format loh_flag
#for input_format, use one of I, S or T. S and T are the preferred options as they work directly from the raw output file from Sequenza or Titan, respectively.
#for the loh_flag, I recommend you set it to 1. This will include all copy-neutral LOH segments and their BAF in the clustering. This can help resolve clonal clusters with few mutations in them. 
#Setting this to zero will ignore LOH events 

args = commandArgs(trailingOnly=TRUE)
print(args)
library(matlab) #this dependency can probably be removed. Note the single line marked that uses a matlab function can probably be replaced with pure R
library(expands)

seg = args[1]

maf = args[2]

sample = args[3]
print(paste("sample:",args[3]))
include_loh = 0

max_score = as.double(args[4])

precision = as.double(args[5])

input_mode = args[6]  # S for sequenza, I for IGV-friendly seg file, T for Titan full segments file, O for augmented OncoSNP .cnvs file

cn_style = args[7] # 1 for integer values, 2 for the actual values based on log ratio (2 is recommended by the EXPANDS author)

# for testing
# seg = "/Volumes/morinlab/projects/2016_dlbcl_dlc_lymphoma_gene_pool/analysis/oncoSNP/3-augmented_oncoSNP/DLC_0010.aug.cnvs"
# maf = "/Volumes/morinlab/projects/2016_dlbcl_dlc_lymphoma_gene_pool/analysis/strelka_pipeline/9-aug_maf/DLC-10_S37.aug.maf"
# #seg = "/Volumes/morinlab/projects/2016_dlbcl_ex_realn/tumour_copy_number/PT003_segments.txt"
# #maf = "/Volumes/morinlab/projects/2016_dlbcl_ex_realn/aug_maf/PT003_Pd15.aug.maf"
# sample = "DLC_0010"
# includ_loh = 1
# max_score = 2.5
# precision = 0.05
# input_mode = "O"
# cn_style = 2

#could also be made an argument if the user wants to really slow down the program by decreasing this :)
min_freq = 0.1
max_PM=5  #6 was causing some weird results

print("expecting input format:")
print(input_mode)

if(length(args)>7){
  #loh flag set to a true value. Default is to not do anything with LOH events
  include_loh = as.double(args[8])
  print(paste("running in LOH mode", include_loh, sep=" "))
  #if set to 1, neutral LOH only will be included, if set to 2, deletion LOH only will be included, if set to 3, all LOH will be included
}
if(input_mode == "T"){
  print(paste("loading Titan output file ",seg))
  #Sample	Chromosome	Start_Position(bp)	End_Position(bp)	Length(bp)	Median_Ratio	Median_logR	TITAN_state	TITAN_call	Copy_Number	MinorCN	MajorCN	Clonal_Cluster	Clonal_Frequency
  seg1=read.csv(seg,stringsAsFactors=FALSE,sep="\t")
  loh_string = "no_LOH_"
  if(include_loh > 0){
    if(include_loh == 1){
      neut = seg1[,"Copy_Number"] == 2
      loh1 = seg1[,"MajorCN"]==2
      keep_loh = neut & loh1
      loh_string = "neut_LOH_"
    } else if(include_loh ==2){
      #warning: These other two options should not be used. I don't think they work well. Please stick with 0 or 1
      del = seg1[,"Copy_Number"] == 1
      loh1 = seg1[,"MajorCN"] + seg1[,"MinorCN"] ==1
      keep_loh = del & loh1
      loh_string = "del_LOH_"
    } else if(include_loh==3){
      loh_string = "any_LOH_"
      not_amp = seg1[,"Copy_Number"] <=2
      loh1 = seg1[,"MinorCN"]==0
      keep_loh = not_amp & loh1
    }
    
    sequenza_loh = seg1[keep_loh,]
    
    vaf_loh_event= sequenza_loh[,"Median_Ratio"]  
    loh_snv_data=matrix(nrow=length(vaf_loh_event),ncol=7,dimnames=list(c(),c("chr","startpos","endpos","REF","ALT","AF_Tumor","PN_B")))
    loh_snv_data[,"chr"] = as.numeric(sequenza_loh[,"Chromosome"])
    loh_snv_data[,"startpos"] = as.numeric(sequenza_loh[,3] + 1)
    loh_snv_data[,"endpos"] = as.numeric(sequenza_loh[,3] + 1)
    loh_snv_data[,"AF_Tumor"]= as.numeric(sequenza_loh[,"Median_Ratio"])
    loh_snv_data[,"REF"] = as.numeric(rep(65,length(vaf_loh_event)))
    loh_snv_data[,"ALT"] = as.numeric(rep(45,length(vaf_loh_event)))
    loh_snv_data[,"PN_B"] = as.numeric(rep(1,length(vaf_loh_event)))
  }
  
  chroms = as.numeric(seg1[,"Chromosome"])
  
  keep_seg = seg1[,4] - seg1[,3] > 1
  
  keep_chrom = !is.na(chroms<23)
  
  keep_both = keep_chrom & keep_seg
  seg2=matrix(nrow=dim(seg1[keep_both,][1]),ncol=4)
  
  colnames(seg2) = c("chr","startpos","endpos","CN_Estimate")
  
  seg2[,"chr"] = as.numeric(seg1[keep_both,"Chromosome"])
  
  seg2[,"startpos"] = as.numeric(seg1[keep_both,3])
  seg2[,"endpos"] = as.numeric(seg1[keep_both,4])
  if(cn_style ==1){
    seg2[,4]=as.numeric(seg1[keep_both,"Copy_Number"])
  }
  else if(cn_style ==2){
    seg2[,4] = 2*2^seg1[keep_both,"Median_logR"]
  }	
  
}else if(input_mode == "S"){
  print(paste("loading Sequenza ",seg))
  seg1=read.csv(seg,stringsAsFactors=FALSE,sep="\t")
  #note that the CNt value is being used directly by Expands. For Titan, one would have to convert the log ratio to absolute copy number
  
  loh_string = "no_LOH_"
  
  if(include_loh > 0){
    if(include_loh == 1){
      neut = seg1[,"CNt"] == 2
      loh1 = seg1[,"A"]==2
      keep_loh = neut & loh1
      loh_string = "neut_LOH_"
    } else if(include_loh ==2){
      #don't use these other two options!
      del = seg1[,"CNt"] == 1
      loh1 = seg1[,"A"] + seg1[,"B"] ==1
      keep_loh = del & loh1
      loh_string = "del_LOH_"
    } else if(include_loh==3){
      loh_string = "any_LOH"
      not_amp = seg1[,"CNt"] <=2
      loh1 = seg1[,"B"]==0
      keep_loh = not_amp & loh1
    }
    
    sequenza_loh = seg1[keep_loh,]
    
    vaf_loh_event= 1- sequenza_loh[,"Bf"]  #the Expands method expects VAFs that were increased relative to the normal but Sequenza reports BAF relative to lower Reference_Allele
    
    loh_snv_data=matrix(nrow=length(vaf_loh_event),ncol=7,dimnames=list(c(),c("chr","startpos","endpos","REF","ALT","AF_Tumor","PN_B")))
    loh_snv_data[,"chr"] = as.numeric(sequenza_loh[,"chromosome"])
    loh_snv_data[,"startpos"] = as.numeric(sequenza_loh[,"start.pos"] + 1)
    loh_snv_data[,"endpos"] = as.numeric(sequenza_loh[,"start.pos"] + 1)
    loh_snv_data[,"AF_Tumor"]= as.numeric(vaf_loh_event)
    loh_snv_data[,"REF"] = as.numeric(rep(65,length(vaf_loh_event)))
    loh_snv_data[,"ALT"] = as.numeric(rep(45,length(vaf_loh_event)))
    loh_snv_data[,"PN_B"] = as.numeric(rep(1,length(vaf_loh_event)))
    #loh_snv_data[,"PN_B"] = as.numeric(rep(0,length(vaf_loh_event))) #don't set as LOH because it works better if you don't. Unsure why
  }
  chroms_a = seg1[,2]
  chroms = as.numeric(chroms_a)
  
  #	chroms = as.numeric(seg1[,"chromosome"])
  
  keep_seg = seg1[,"end.pos"] - seg1[,"start.pos"] > 1
  
  keep_chrom = !is.na(chroms<23)
  
  keep_both = keep_chrom & keep_seg
  seg2=matrix(nrow=dim(seg1[keep_both,][1]),ncol=4)
  
  colnames(seg2) = c("chr","startpos","endpos","CN_Estimate")
  
  seg2[,"chr"] = as.numeric(seg1[keep_both,"chromosome"])
  
  seg2[,"startpos"] = as.numeric(seg1[keep_both,"start.pos"])
  seg2[,"endpos"] = as.numeric(seg1[keep_both,"end.pos"])
  if(cn_style ==1){
    seg2[,4]=as.numeric(seg1[keep_both,"CNt"])
  } else if (cn_style == 2){
    logratio=log2(as.numeric(seg1[keep_both,"depth.ratio"]))
    seg2[,4] = 2*2^logratio
  }
  
} else if (input_mode == "I"){
  #LOH information needs to be included in these files and BAF
  #example format: sample	chr	start	end	LOH_flag	BAF	median
  seg1=read.csv(seg,stringsAsFactors=FALSE,sep="\t")
  #convert last column from log ratio to absolute CN
  logratios = seg1[,length(seg1[1,])]
  abs = 2*2^logratios
  chroms_a = seg1[,2]
  chroms = as.numeric(chroms_a)
  keep_chrom= !is.na(chroms<23)
  
  #loh_string = "no_LOH_"
  
  loh_status = seg1[,"LOH_flag"]
  
  
  if(include_loh > 0){
    #Warning: This is the most experimental mode of running currently. 
    
    sequenza_loh = seg1[loh_status ==1,]
    vaf_loh_event= 1- sequenza_loh[,"BAF"]  #the Expands method expects VAFs that were increased relative to the normal but Sequenza reports BAF relative to lower Reference_Allel
    
    loh_snv_data=matrix(nrow=length(vaf_loh_event),ncol=7,dimnames=list(c(),c("chr","startpos","endpos","REF","ALT","AF_Tumor","PN_B")))
    loh_snv_data[,"chr"] = as.numeric(sequenza_loh[,"chromosome"])
    loh_snv_data[,"startpos"] = as.numeric(sequenza_loh[,"start"] + 1)
    loh_snv_data[,"endpos"] = as.numeric(sequenza_loh[,"end"] + 1)
    loh_snv_data[,"AF_Tumor"]= as.numeric(vaf_loh_event)
    loh_snv_data[,"REF"] = as.numeric(rep(65,length(vaf_loh_event)))
    loh_snv_data[,"ALT"] = as.numeric(rep(45,length(vaf_loh_event)))
    loh_snv_data[,"PN_B"] = as.numeric(rep(1,length(vaf_loh_event)))
    #print(loh_snv_data)
    #q()
  }
  seg2=matrix(nrow=length(seg1[keep_chrom,1]),ncol=4)
  
  #seg2=matrix(nrow=length(seg1[keep_chrom,1]),ncol=4)
  colnames(seg2) = c("chr","startpos","endpos","CN_Estimate")
  
  seg2[,"chr"] = as.numeric(seg1[keep_chrom,2])
  
  seg2[,"startpos"] = as.numeric(seg1[keep_chrom,3])
  seg2[,"endpos"] = as.numeric(seg1[keep_chrom,4])
  seg2[,4]=abs[keep_chrom]
  
} else if (input_mode == "O") {
  
  # Expects OncoSNP .cnvs file augmented with Log R Ratio and BAF using oncosnputils
  # Columns: chr	start	end	copyNum	loh	rank	logLik	numMarkers	normFrac	state	ploidyNum
  # majorCopyNumber minorCopyNumber	LRR	LRRShifted	BAF	numProbes	numSnpProbes
  
  print(paste("loading OncoSNP ", seg))
  seg1 <- read.csv(seg, stringsAsFactors = FALSE, sep = "\t")
  chroms_a  <- seg1[, 1]
  chroms <- as.numeric(chroms_a)
  
  keep_seg <- seg1[, "end"] - seg1[, "start"] > 1
  keep_chrom  <- !is.na(chroms < 23)
  
  keep_both <- keep_chrom & keep_seg
  seg2 <- matrix(nrow = dim(seg1[keep_both, ][1]), ncol = 4)
  
  colnames(seg2) <- c("chr", "startpos", "endpos", "CN_Estimate")
  
  seg2[, "chr"] <- as.numeric(seg1[keep_both, "chr"])
  seg2[, "startpos"] <- as.numeric(seg1[keep_both, "start"])
  seg2[, "endpos"] <- as.numeric(seg1[keep_both, "end"])
  
  if (cn_style == 1) {
    
    seg2[, 4] <- as.numeric(seg1[keep_both, "copyNum"])
    
  } else if (cn_style == 2) {
    
    logratio <- as.numeric(seg1[keep_both, "LRR"])
    seg2[, 4] <- 2*2^logratio
    
  }
  
  
} else {
  print("no input mode specified!")
  q()
}


#this code will likely be removed as the new release seems to work better for deletions

mask_deletions = 0 

if(cn_style == 2){
  mask_deletions = 0
}

if(mask_deletions){
  seg2.mask = seg2[,4] < 2
  seg2.save = seg2[,"CN_Estimate"]
  
  seg2[seg2.mask,"CN_Estimate"] = 2
  seg2=cbind(seg2,seg2.save)
  colnames(seg2)= c("chr","startpos","endpos","CN_Estimate","CN_Estimate_nomask")
  
}

print(paste("loading MAF ",maf))
maf_data = read.csv(maf,sep="\t",stringsAsFactors=FALSE )


maf_keep_cols = c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","t_ref_count","t_alt_count","n_ref_count","n_alt_count")

#this method works with Indels but is a bit hacky because a fake ref/alt value is used. It doesn't really matter as long as you ignore those values in the output. Expands ignores them AFAIK
maf_keep = maf_data[,maf_keep_cols]
maf_keep[,"Reference_Allele"]=sapply(maf_keep[,"Reference_Allele"],function(x) substr(x,1,1))
maf_keep[,"Tumor_Seq_Allele2"]=sapply(maf_keep[,"Tumor_Seq_Allele2"],function(x) substr(x,1,1))

snv_data=matrix(nrow=dim(maf_keep[1]),ncol=7,dimnames=list(c(),c("chr","startpos","endpos","REF","ALT","AF_Tumor","PN_B")))

# Remove chr prefix
snv_data[,"chr"] <- sub("^chr", "", maf_data[,"Chromosome"])


#load start and end position into matrix
snv_data[,"startpos"] = as.numeric(maf_keep[,"Start_Position"])
snv_data[,"endpos"] = as.numeric(maf_keep[,"End_Position"])
snv_data[,"REF"]=sapply(maf_keep[,"Reference_Allele"],utf8ToInt)
snv_data[,"ALT"]=sapply(maf_keep[,"Tumor_Seq_Allele2"],utf8ToInt)
snv_data[,"AF_Tumor"]=maf_keep[,"t_alt_count"]/(maf_keep[,"t_alt_count"]+maf_keep[,"t_ref_count"])

#set flag to somatic for all SNVs
snv_data[,"PN_B"] = 0


#this code can all be removed for the Galaxy tool version ending with <-HERE
#for convenience, make a pyclone input file using some code recycled from expands. There might be a more reasonable place to declare this function
assignStatesToMutation<-function(dm,cbs,cols){
  print("Assigning copy number to mutations for PyClone...")
  ##Assign copy numbers in cbs to mutations in dm
  for (k in 1:nrow(cbs)){
    
    idx=which(dm[,"chr"]==cbs[k,"chr"] & dm[,"startpos"]>=cbs[k,"startpos"] & dm[,"startpos"]<=cbs[k,"endpos"]);
    if (length(idx)==0){
      next;
    }
    #the next line is the only bit of code that relies on the matlab dependencie
    dm[idx,cols]=repmat(as.numeric(cbs[k,cols]),length(idx),1);
    dm[idx,"normal_cn"] = rep(2,length(idx))
  }
  dm=dm[,colnames(dm)!="segmentLength"];
  print("... Done.")
  
  return(dm);
}


py_snv_data=matrix(nrow=dim(maf_keep[1]),ncol=10,dimnames=list(c(),c("gene","chr","startpos","endpos","mutation_id","ref_counts","var_counts","normal_cn","major_cn","minor_cn")))

# Remove chr prefix
py_snv_data[, "chr"] <- sub("^chr", "", maf_keep[,"Chromosome"])

#load start and end position into matrix
py_snv_data[,"startpos"] = as.numeric(maf_keep[,"Start_Position"])
py_snv_data[,"endpos"] = as.numeric(maf_keep[,"End_Position"])

py_snv_data[,"ref_counts"]=maf_keep[,"t_ref_count"]
py_snv_data[,"var_counts"]=maf_keep[,"t_alt_count"]

if(input_mode == "S") {
  
  #rename columns in Sequenza data to match normal_cn, minor_cn, major_cn and position/chromoosome name style of Expands
  colnames(seg1)=c("chr","startpos","endpos","Bf","N.BAF","sd.BAF","depth.ratio","N.ratio","sd.ratio","CNt","major_cn","minor_cn","segmentLength") #note, we need to fill normal_cn with 2 for everything
  #set segment length
  seg1[,"segmentLength"] = seg1[,"endpos"] - seg1[,"startpos"]
  py_snv_data_assigned=assignStatesToMutation(py_snv_data,seg1,c("minor_cn","major_cn"))
  py_snv_data_assigned[,"gene"]=maf_keep[,"Hugo_Symbol"]
  py_snv_data_assigned[,"mutation_id"]=paste(py_snv_data_assigned[,"gene"],py_snv_data_assigned[,"startpos"],sep="_")
  out_pyclone = paste("./",sample,"_pyclone_in.tsv",sep="")
  keepers = !is.na(py_snv_data_assigned[,"normal_cn"])
  write.table(py_snv_data_assigned[keepers,],file=out_pyclone,sep="\t",quote=FALSE)
  
} else if (input_mode == "O") {
  
  colnames(seg1) <- c("chr", "startpos", "endpos", "CN", "loh", "rank", "logLik", "numMarkers", 
                     "normFrac", "state", "ploidyNum", "major_cn", "minor_cn", "log.ratio", "log.ratio.shifted", "BAF", "numProbes", "numSnpProbes")
  seg1[, "normal_cn"] <- 2
  seg1[, "segmentLength"] <- seg1[, "endpos"] - seg1[, "startpos"]
  py_snv_data_assigned <- assignStatesToMutation(py_snv_data, seg1, c("minor_cn","major_cn"))
  py_snv_data_assigned[, "gene"] <- maf_keep[, "Hugo_Symbol"]
  py_snv_data_assigned[, "mutation_id"] <- paste(py_snv_data_assigned[, "gene"],
                                                py_snv_data_assigned[, "startpos"], sep = "_")
  
  out_pyclone <- paste0("./", sample, "_pyclone_in.tsv")
  keepers <- !is.na(py_snv_data_assigned[, "normal_cn"])
  write.table(py_snv_data_assigned[keepers, ], file = out_pyclone, sep = "\t", quote = FALSE)

}



#<-HERE

if(include_loh > 0){
  
  merge_snv=rbind(loh_snv_data,snv_data)
  #loh_string = "LOH_"
  samp_param = paste(sample,loh_string,"_state_INDEL_DelMask_maxpm_", max_PM, "_score_", max_score, ",precision_", precision,sep="")
} else{
  loh_string = "no_LOH"
  merge_snv = snv_data
  samp_param = paste(sample,loh_string,"_noLOH_INDEL_DelMask_maxpm_", max_PM, "_score_", max_score, ",precision_", precision,sep="")
}
print("merged")
print(merge_snv)
print(".....")
print(seg2)


#run individual expands steps separately
dm=assignQuantityToMutation(merge_snv,seg2,"CN_Estimate")

plotF=1

ii=which(is.na(dm[,"CN_Estimate"]));
if (length(ii)>0){
  print(paste(length(ii), " SNV(s) excluded due to unavailable copy number in that region."));
  dm=dm[-ii,];
}
ii=which(dm[,"CN_Estimate"]<1);
homDelRegions=c();
if (length(ii)>0){
  print(paste(length(ii), " SNV(s) excluded due to homozygous deletions within that region."));
  homDelRegions=dm[ii,];
  dm=dm[-ii,];
}
ii=which(dm[,"CN_Estimate"]>max_PM);
if (length(ii)>0){
  print(paste(length(ii), " SNV(s) excluded due to high-level amplifications (>",max_PM, "copies) within that region. Consider increasing value of parameter max_PM to facilitate inclusion of these SNVs, provided high coverage data (> 150 fold) is available"));
  dm=dm[-ii,];
}
ii=which(dm[,"AF_Tumor"]*dm[,"CN_Estimate"]<min_freq);
if (length(ii)>0){
  print(paste(length(ii), " SNV(s) excluded due to AF*CN below ", min_freq," (SNV can't be explained by an SP present in ",min_freq*100 ,"% or more of the sample)."));
  dm=dm[-ii,];
}

if(precision ==1){
  precision=0.1/log(nrow(dm)/7); #use default?
  
}
print(dm)
print(loh_snv_data)
print(samp_param)
dirF=getwd();

output_file=paste(dirF, "/", samp_param,"unmodified.sps",sep="");

file = paste(samp_param,".pdf",sep="")

cfd=computeCellFrequencyDistributions(dm, max_PM, precision, min_CellFreq=min_freq)

toUseIdx=which(apply(is.finite(cfd$densities),1,all) )


SPs=clusterCellFrequencies(cfd$densities[toUseIdx,], precision, min_CellFreq=min_freq)


SPs=SPs[SPs[,"score"]<=max_score,] 

print("subpopulation details:")

print(SPs)

print("running assignMutations...")



aM= assignMutations( dm, SPs,max_PM=max_PM)


#unmaks the deletions for vizualization
#can probably be removed
if(mask_deletions){
  seg2[,"CN_Estimate"] = seg2[,"CN_Estimate_nomask"]
  dm = assignQuantityToMutation(merge_snv,seg2,"CN_Estimate")
  some=dm[,"startpos"] %in% aM$dm[,"startpos"]
  
  aM$dm[,"CN_Estimate"] = dm[some,"CN_Estimate"]
  
}
#save SPS file from Expands with mutations assigned to subclones


#file = "expands_plot_simu_incl_chr6LOH.pdf"
#plot the Expands image showing mutations assigned to their SPs

#aM$dm[,"%maxP"] = 1

file1 = paste("rawPlot_",file,sep="")
pdf(file1)
plotSPs(aM$dm,sampleID=sample,cex=1,rawAF=TRUE)
dev.off()

pdf(file)
plotSPs(aM$dm, sampleID=sample,cex=1)
dev.off()

#print(dm)

print("final SPs")
print(aM$finalSPs)

#we should probably write out the SP details to another file as is done in the runExpands wrapper

suppressWarnings(write.table(aM$dm,file = output_file, quote = FALSE, sep = "\t", row.names=FALSE))

