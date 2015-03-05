#BEDTOOLS FUNCTION DEFINITIONS
BED.COLS <- c("chr","start","end")
#~~~~# Code from Rosemary McCloskey
## FUNCTION DEFINITIONS FOR BED FILE ANALYSIS 
# write a data frame to a BED file
# data frame must have the columns "chr", "start", and "end"
# return the file name written to
write.bed <- function (bed.df, bed.file=tempfile()) {
  stopifnot(colnames(bed.df)[1:length(BED.COLS)] == BED.COLS)
  write.table(bed.df, bed.file, col.names=F, row.name=F, sep="\t", quote=F)
  bed.file
}

# run a bedtools command 
# bed files should both have chr, start, end as their first three columns
bedtools <- function (command, bed.a, bed.b, args="") {
  bed.file.a <- write.bed(bed.a)
  cmd <- paste("bedtools", command, args)
  if (command == "sort") {
    bed.file.a <- write.bed(bed.a)
    cmd <- paste(cmd,"-i", bed.file.a)
  } else{
  if (!is.null(bed.b)) {
    bed.file.b <- write.bed(bed.b)
    cmd <- paste(cmd, "-a", bed.file.a, "-b", bed.file.b)
  } else {
    cmd <- paste(cmd, "-i", bed.file.a)
  }
  }
  cat("Running", cmd, "... ")
  output <- paste(system(cmd, intern=T), collapse="\n")
  cat("done\n")
  if (output == "") return (NULL)
  res <- read.table(textConnection(output))
  if (ncol(res) == 3)
    setNames(res, BED.COLS)
  else if (ncol(res) == ncol(bed.a))
    setNames(res, colnames(bed.a))
  else if (ncol(res) == ncol(bed.a) + ncol(bed.b))
    setNames(res, make.names(c(colnames(bed.a), colnames(bed.b)), unique=T))
  else if (ncol(res) == ncol(bed.a) + ncol(bed.b) + 1)
    setNames(res, make.names(c(colnames(bed.a), colnames(bed.b), "overlap"), unique=T))
  else
    res 
}

#Reformat sample id's to match the format defined in data files
string_trimmer <- function (input_name) {
  filename=as.character(input_name)
  filename = gsub("(^[[:space:]]+|[[:space:]]+$)", "", filename)  
  if((substring(filename,nchar(filename)-1,nchar(filename)) > 18) & ((substring(filename,1,2) == "PT")==TRUE)){
    filename=paste(substring(filename,1,2),substring(filename,nchar(filename)-1,nchar(filename)),sep="")
  }
  input_name = filename
  input_name
}

#Get CCF using variant allele fraction, purity, copy number, prevalence
vaf.to.ccf <- function (vaf, pur, cn, prev) {
  avg.copies <- (cn*prev + 2*(1-prev))*pur + 2*(1-pur)
  pur.bound <- (cn-1)*prev*pur/avg.copies
  ccf.bound <- pur/avg.copies
  midpoint <- mean(c(pur.bound, ccf.bound))
  ifelse(cn <= 2 | vaf <= midpoint,
         avg.copies*vaf/pur,
         avg.copies*vaf/pur - (cn-2)*prev)
}

#~~~~# End of code from Rosemary McCloskey

#DIRECTORIES
homedir = "/Volumes/Shared/projects/" #Location hosting dlbcl project
home_titan = "dlbcl_relapse_exome/seg_files/TITAN/" #Location in dlbcl project, hosting TITAN output data
snp_txt = "/Users/jgrewal/Desktop/CCF/snp.txt" #Variant input file in the format chr start stop patient
localdir = "/Users/jgrewal/Desktop/CCF/"
snp_complete="/Users/jgrewal/Desktop/CCF/snp_complete.txt"
#read in SNP variant file and extract patient ids
snp = read.delim(snp_txt,header=FALSE,stringsAsFactors=FALSE) #patient, chr, startpos, endpos, x, x, median_log, x, x, copy_number
#sample_id = unique(snp[4]) #list of unique sample ids
filenames = snp[4]
for (i in 1:dim(filenames)[1]){
  filenames[i,1] = string_trimmer(filenames[i,1])
}
sample_id = unique(filenames)
rownames(sample_id) <- NULL
snp['samples'] = snp[4]
snp[4] = filenames

#initialize holder dataframe
bedresult_cols = c("chr_var","start_var","end_var","chr_cnv","start_cnv","end_cnv","overlap","TITAN_call","local_cnv","prevalence","purity","patient")
bedresult =  data.frame(matrix(vector(), 0, length(bedresult_cols), dimnames=list(c(), bedresult_cols)), stringsAsFactors=F)

#Process each patient
for (j in 1:dim(sample_id)[1]){
  PROCESSFLAG = 1
  patient = as.character(sample_id[j,1])
  #get correct file
  filename = as.character(patient)
  pt_titan = paste(homedir,home_titan,filename,"/",sep="")
  
  #If params file exists, get the purity estimate
  paramfile = list.files(path=pt_titan, pattern= "*_params.txt")
  path_paramfile=paste(pt_titan,paramfile,sep="")
  if(file.exists(path_paramfile)){
    temp = read.delim(path_paramfile,nrow=1,header=FALSE,stringsAsFactors=FALSE)
    purity = 1-temp[1,2]
  }else{
    PROCESSFLAG = 0
    print (paste("Params file for patient ",patient," not found",sep=""))
  }
  
  #If Titan's cnv output file exists, get the segments of CN variation from the file
  #Also get the prevalence values (for CCF calculation)
  cnvfile = list.files(path=pt_titan, pattern= "*_segs.txt")
  path_cnvfile=paste(pt_titan,cnvfile,sep="")
  if(file.exists(path_cnvfile)){
    patient_cnv = read.table(pipe(paste("cut -f 2,3,4,10",path_cnvfile,sep=" ")),header=TRUE,stringsAsFactors=FALSE) #chr start end copy_number
    rownames(patient_cnv) <- NULL
    colnames(patient_cnv) <- c("chr","start","end","local_cnv")
    patient_prevalence = read.table(pipe(paste("cut -f 2,3,4,10,14,9",path_cnvfile,sep=" ")),header=TRUE,stringsAsFactors=FALSE) #chr start end copy_number cellular frequency (i.e. prevalence)
  }else{
    PROCESSFLAG = 0
    print (paste("Titan seg file for patient ",patient," not found",sep=""))
  }
  
  #Create patient specific snp bed data
  patient_snp = snp[snp[4] == patient,1:3]
  rownames(patient_snp) <- NULL
  colnames(patient_snp) <- c("chr","start","end")
  
  #Prepare output bed files (patient_snp.bed, patient_cnv_exons.bed)
  if((dim(patient_snp)[1]>0) & PROCESSFLAG){
    temp2 <- transform(patient_snp,chr=as.integer(chr))
    temp_snp_table = temp2[order(temp2[1],decreasing=FALSE),]
    rownames(temp_snp_table) <- NULL
    temp_snp_table = bedtools("sort",patient_snp)
    temp_cnv_table = bedtools("sort",patient_cnv[1:3])
    #temp_table<-bedtools("closest",temp_snp_table,patient_cnv[1:3],"-d")
    temp_table<-bedtools("closest",temp_snp_table,temp_cnv_table,"-d")
    temp_table<-temp_table[temp_table['overlap']>-1,1:7] #Select segments with overlap
    rownames(temp_table) <- NULL
    #get copy numbers in these regions, that have mapped close/on to variant position
    for (i in 1:dim(temp_table)[1]){
      temp_table[i,"TITAN_call"] =  patient_prevalence[do.call(paste, patient_prevalence[1:3]) %in% paste(c(paste(temp_table[i,'chr.1']), temp_table[i,'start.1'], temp_table[i,'end.1']), collapse = " "),"TITAN_call"]
      temp_table[i,"local_cnv"] = patient_cnv[do.call(paste, patient_cnv[1:3]) %in% paste(c(paste(temp_table[i,'chr.1']), temp_table[i,'start.1'], temp_table[i,'end.1']), collapse = " "),4]
      temp_table[i,"prevalence"] = patient_prevalence[do.call(paste, patient_prevalence[1:3]) %in% paste(c(paste(temp_table[i,'chr.1']), temp_table[i,'start.1'], temp_table[i,'end.1']), collapse = " "),"Clonal_Frequency"]
    }
    temp_table[,"purity"] = purity
    temp_table[,"patient"] = patient
    colnames(temp_table) <- bedresult_cols
    bedresult <- rbind(bedresult, temp_table)
  }
}

#Print resulting file with snv locations and matching cnv segments
bedresult<-unique(bedresult)
write.table(bedresult,paste(localdir,"local_cn_result.txt",sep=""),sep="\t",row.names=F)
snp_whole = read.delim(snp_complete,header=TRUE,stringsAsFactors=FALSE) #patient, chr, startpos, endpos, x, x, median_log, x, x, copy_number

for (i in 1:dim(snp_whole)[1]){
  patient_temp = snp_whole[i,'SampleID']
  patient_temp = string_trimmer(patient_temp)
  #cat('(', i,')',patient_temp,'..')
  #If patient's bed matches exist
  if(dim(bedresult[bedresult['patient']==patient_temp,])[1] > 0){
  #Append purity to corresponding patient variant record
  snp_whole[i,"Purity.estimate"] = unique(bedresult[bedresult['patient']==patient_temp,"purity"])
  #Append local tital call to corresponding patient variant record
  temp_call = bedresult[do.call(paste,bedresult[c("chr_var","start_var","end_var","patient")]) %in% paste(c(paste(snp_whole[i,'Chromosome']), snp_whole[i,'Start_Position'], snp_whole[i,'End_Position'],patient_temp), collapse = " "),"TITAN_call"]
  snp_whole[i,"TITAN_call"] = ifelse(length(temp_call)>0, temp_call, NA)
  #Append local cnv estimate to corresponding patient variant record
  temp_cnv = bedresult[do.call(paste,bedresult[c("chr_var","start_var","end_var","patient")]) %in% paste(c(paste(snp_whole[i,'Chromosome']), snp_whole[i,'Start_Position'], snp_whole[i,'End_Position'],patient_temp), collapse = " "),"local_cnv"]
  snp_whole[i,"Local.copy.number"] = ifelse(length(temp_cnv)>0, temp_cnv, NA)
  #Append local prevalence estimate to corresponding patient variant record
  temp_prev = bedresult[do.call(paste,bedresult[c("chr_var","start_var","end_var","patient")]) %in% paste(c(paste(snp_whole[i,'Chromosome']), snp_whole[i,'Start_Position'], snp_whole[i,'End_Position'],patient_temp), collapse = " "),"prevalence"]
  snp_whole[i,"Prevalence.estimate"] = ifelse(length(temp_prev)>0, temp_prev, NA)
  #Append CCF calculation result to corresponding patient variant record
  allele_frac = as.numeric(snp_whole[i,"Variant.Allele.fraction"])
  purity = snp_whole[i,"Purity.estimate"]
  copy_num = snp_whole[i,"Local.copy.number"]
  preval = 1 #ifelse(is.na(snp_whole[i,"Prevalence.estimate"]) & (snp_whole[i,"TITAN_call"]=="HET"),1,snp_whole[i,"Prevalence.estimate"])
  snp_whole[i,"Prevalence.estimate"] = preval
  snp_whole[i,"CCF.estimate"] = ifelse(is.na(preval),NA,vaf.to.ccf(allele_frac,purity,copy_num,preval))
  }
}

write.table(snp_whole,paste(localdir,"resultfile_snp.txt",sep=""),sep="\t",row.names=F)


#Post preivew
as.data.frame(table(snp[4]))[order(as.data.frame(table(snp[4]))[1]),]
as.data.frame(table(bedresult['patient']))[order(as.data.frame(table(bedresult['patient']))[1]),]



