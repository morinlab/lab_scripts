### USAGE ##########################################################################################
usage <- function() {
  
  usage.text <-
    '\nNAME
  ccf_loh.R is an R script for calculating ccf for multiple samples.
  
  USAGE
  Rscript splitBAM.R [OPTIONS]
  
  OPTIONS
  -h, --help
  -s, --snvfile (optional)
  The consolidated snv text file for multiple samples* with the correct header**. 
    [Default is /Volumes/morinlab/projects/dlbcl_relapse_exome/ccf/data/myinfile_snp.txt]
  -b, --output_bed (required)
  The destination for the bed file. 
    [Default is ./tempbedfile.txt]
  -o, --output_snv (required)
  The destination for the final results file. 
    [Default is ./snv_results.txt]
  -t, --titandir (optional)
  The default directory for Titan results, organized*** by sample* 
    [Default is /Volumes/morinlab/projects/dlbcl_relapse_exome/seg_files/TITAN/]
  -v, --vaf_mean_file (optional)
  Txt file containing sampleid* and meanVAF, to be used to estimate purity.
    [Default is /Volumes/morinlab/projects/dlbcl_relapse_exome/ccf/data/vaf_dlbcl.txt]

  INPUT FORMAT
  *   The sample names must be consistent between the input snv file, the Titan per-sample directories, and the meanVAF text file.
  **  The input snv file must have the following 6 column names in the header
        - \'SampleID\'
        - \'Chromosome\'
        - \'Start_Position\'
        - \'End_Position\'
        - \'Variant Allele fraction\'
        - \'Purity estimate\'
  *** The directory structure for your TITAN results must abide by the following format:
        - SampleID1*
            ->contains 1 _params.txt and 1 _segs.txt file (based on optimal cluster)
        - SampleID2*
            ->contains 1 _params.txt and 1 _segs.txt file (based on optimal cluster)
        ...and so on for subsequent samples
  EXAMPLE USAGE
  Rscript  splitBAM.R -s mysnvfile.txt -b mybedfile.bed -o mysnvfile_withccf.txt
  
  AUTHOR
  Jasleen Grewal\n\n';

	return(usage.text);
	}

#DIRECTORIES
#temp <- readLines("/Volumes/morinlab/projects/dlbcl_relapse_exome/ccf/directories.txt")
#eval(parse(text=temp))
suppressWarnings(suppressMessages(library(getopt,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))
spec=matrix(c(
  'help',  		'h', 0, "logical",
  'snvfile','s',0,"character", #"The consolidated snv txt file.",
  'output_bed','b',1,"character", # "Destination for bed file.",
  'output_snv','o',1,"character", # "Destination for snv file with CCF results. Final results file.",
  'titandir','t',0,"character", # "Path from home directory to titan results",
  'vaf_mean_file',"v",0,"character" # "Txt file containing sampleid and meanVAF, to be used to estimate purity."
  ), byrow=TRUE,ncol=4)
opt=getopt(spec);

#Set defaults
if (is.null(opt$snvfile)) {opt$snvfile = "/Volumes/morinlab/projects/dlbcl_relapse_exome/ccf/data/myinfile_snp.txt"}
if (is.null(opt$vaf_mean_file)) {opt$vaf_mean_file = "/Volumes/morinlab/projects/dlbcl_relapse_exome/ccf/data/vaf_dlbcl.txt"}
if (is.null(opt$titandir)) {opt$titandir = "/Volumes/morinlab/projects/dlbcl_relapse_exome/seg_files/TITAN/"}
if (is.null(opt$output_bed)) {opt$output_bed = "tempbedfile.bed"}
if (is.null(opt$output_snv)) {opt$output_snv = "snv_results.txt"}

if( !is.null(opt[['help']]) || is.null(opt[['snvfile']])) {
  cat(usage());
  q(status=1);
}

home_titan=opt$titandir
snp_complete_txt=opt$snvfile
vaf_recalculated_txt=opt$vaf_mean_file
bed_output_loh=opt$output_bed
results_output_loh=opt$output_snv

#----------------------------------------------------------------------------------#

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

#Get CCF using variant allele fraction, purity, copy number, prevalence
vaf.to.ccf <- function (vaf, pur, prev, minor.cn, major.cn) {
  cn <- major.cn + minor.cn
  alpha <- (cn*prev + 2*(1-prev))*pur + 2*(1-pur)
  # no duplications
  if (minor.cn <= 1 & major.cn <= 1) {
    alpha*vaf/pur
    # one duplication, no LOH
  } else if (minor.cn == 1) {
    if (vaf >= (major.cn*prev+1)*pur/(2*alpha))
      alpha*vaf/pur - (major.cn-1)*prev
    else
      alpha*vaf/pur
    # one duplication with LOH
  } else if (minor.cn == 0) {
    if (vaf >= (1+(major.cn-1)*prev)*pur/(2*alpha))
      alpha*vaf/pur - (major.cn-1)*prev
    else
      alpha*vaf/pur
    # two duplications
  } else {
    if (vaf <= (1+(minor.cn-1)*prev)*pur/(2*alpha))
      alpha*vaf/pur
    else if (vaf >= (1+(major.cn+minor.cn-1)*prev)*pur/(2*alpha))
      alpha*vaf/pur - (major.cn-1)*prev
    else
      alpha*vaf/pur - (minor.cn-1)*prev
  } 
}

#~~~~# End of code from Rosemary McCloskey

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

#Read data
snp_complete = read.delim(snp_complete_txt,header=TRUE,stringsAsFactors=FALSE)
vaf_recalculated = read.delim(vaf_recalculated_txt,header=TRUE,stringsAsFactors=FALSE)
#read in SNP variant file and extract patient ids
filenames = snp_complete['SampleID']
for (i in 1:dim(filenames)[1]){
  filenames[i,1] = string_trimmer(filenames[i,1])
}
sample_id = unique(filenames)
rownames(sample_id) <- NULL
#snp['samples'] = snp[4]
#snp[4] = filenames
#snp_complete['SampleID'] = filenames
#initialize holder dataframe
bedresult_cols = c("chr_var","start_var","end_var","chr_cnv","start_cnv","end_cnv","overlap","TITAN_call","local_cnv","minor_cn","major_cn","prevalence","input_purity","titan_purity","patient","purity_from_mean","vaf_mean")
bedresult =  data.frame(matrix(vector(), 0, length(bedresult_cols), dimnames=list(c(), bedresult_cols)), stringsAsFactors=F)

#Process each patient
for (j in 1:dim(sample_id)[1]){
  PROCESSFLAG = 1
  patient = as.character(sample_id[j,1])
  #get correct file
  filename = as.character(patient)
  pt_titan = paste(home_titan,filename,"/",sep="")
  
  #If params file exists, get the purity estimate
  paramfile = list.files(path=pt_titan, pattern= "*_params.txt")
  path_paramfile=paste(pt_titan,paramfile,sep="")
  if(file.exists(path_paramfile)){
    temp = read.delim(path_paramfile,nrow=1,header=FALSE,stringsAsFactors=FALSE)
    purity = 1-temp[1,2]
  }else{
    PROCESSFLAG = 0
    purity = NA
    print (paste("Params file for patient ",patient," not found",sep=""))
  }
  
  #If Titan's cnv output file exists, get the segments of CN variation from the file
  #Also get the prevalence values, major and minor CN (for CCF calculation)
  cnvfile = list.files(path=pt_titan, pattern= "*_segs.txt")
  path_cnvfile=paste(pt_titan,cnvfile,sep="")
  if(file.exists(path_cnvfile)){
    patient_cnv = read.table(pipe(paste("cut -f 2,3,4,10,11,12",path_cnvfile,sep=" ")),header=TRUE,stringsAsFactors=FALSE) #chr start end copy_number
    rownames(patient_cnv) <- NULL
    colnames(patient_cnv) <- c("chr","start","end","local_cnv","minor_cn","major_cn")
    patient_prevalence = read.table(pipe(paste("cut -f 2,3,4,10,14,9",path_cnvfile,sep=" ")),header=TRUE,stringsAsFactors=FALSE) #chr start end copy_number cellular frequency (i.e. prevalence)
  }else{
    PROCESSFLAG = 0
    print (paste("Titan seg file for patient ",patient," not found",sep=""))
  }
  
  #Create patient specific snp bed data
  #######*SNP DATA FILE (maf)*##############
  provided_purity= unique( snp_complete[snp_complete['SampleID'] == filename,"Purity.estimate"])
  patient_snp = snp_complete[snp_complete['SampleID'] == filename,c("Chromosome","Start_Position","End_Position")] #assume snp data file's first three cols are chr start end
  rownames(patient_snp) <- NULL
  colnames(patient_snp) <- c("chr","start","end")
  
  #Prepare output bed files (patient_snp.bed, patient_cnv_exons.bed)
  if((dim(patient_snp)[1]>0) & PROCESSFLAG){
    temp_snp_table = bedtools("sort",patient_snp)
    temp_cnv_table = bedtools("sort",patient_cnv[1:3])
    #temp_table<-bedtools("closest",temp_snp_table,patient_cnv[1:3],"-d")
    temp_table<-bedtools("closest",temp_snp_table,temp_cnv_table,"-d")
    temp_table<-temp_table[temp_table['overlap']>-1,1:7] #Select segments with overlap
    rownames(temp_table) <- NULL
    #get copy numbers in these regions, that have mapped close/on to variant position
    for (i in 1:dim(temp_table)[1]){
      temp_table[i,"TITAN_call"] =  patient_prevalence[do.call(paste, patient_prevalence[1:3]) %in% paste(c(paste(temp_table[i,'chr.1']), temp_table[i,'start.1'], temp_table[i,'end.1']), collapse = " "),"TITAN_call"]
      temp_table[i,"local_cnv"] = patient_cnv[do.call(paste, patient_cnv[1:3]) %in% paste(c(paste(temp_table[i,'chr.1']), temp_table[i,'start.1'], temp_table[i,'end.1']), collapse = " "),"local_cnv"]
      temp_table[i,"minor_cn"] = patient_cnv[do.call(paste, patient_cnv[1:3]) %in% paste(c(paste(temp_table[i,'chr.1']), temp_table[i,'start.1'], temp_table[i,'end.1']), collapse = " "),"minor_cn"]
      temp_table[i,"major_cn"] = patient_cnv[do.call(paste, patient_cnv[1:3]) %in% paste(c(paste(temp_table[i,'chr.1']), temp_table[i,'start.1'], temp_table[i,'end.1']), collapse = " "),"major_cn"]
      temp_table[i,"prevalence"] = patient_prevalence[do.call(paste, patient_prevalence[1:3]) %in% paste(c(paste(temp_table[i,'chr.1']), temp_table[i,'start.1'], temp_table[i,'end.1']), collapse = " "),"Clonal_Frequency"]
    }
    temp_table[,"input_purity"] = provided_purity
    temp_table[,"titan_purity"] = purity
    temp_table[,"patient"] = patient
    temp_table[,"purity_from_mean"] =(mean(as.numeric(snp_complete[snp_complete["SampleID"] == patient, "Variant.Allele.fraction"])))/0.5
    temp_table[,"vaf_mean"] = as.numeric(vaf_recalculated[vaf_recalculated['patient'] == patient,'vaf'])
    colnames(temp_table) <- bedresult_cols
    bedresult <- rbind(bedresult, temp_table)
  }
}

#Print resulting file with snv locations and matching cnv segments
bedresult<-unique(bedresult)
write.table(bedresult,bed_output_loh,sep="\t",row.names=F)

#CALCULATE CCFs
snp_whole = read.delim(snp_complete_txt,header=TRUE,stringsAsFactors=FALSE) #patient, chr, startpos, endpos, x, x, median_log, x, x, copy_number
#snp_whole = within(snp_complete,rm("CCF.estimate"))
snp_whole=snp_complete
for (i in 1:dim(snp_whole)[1]){
  patient_temp = snp_whole[i,'SampleID']
  patient_temp = string_trimmer(patient_temp)
  #cat('(', i,')',patient_temp,'..')
  #If patient's bed matches exist
  if(dim(bedresult[bedresult['patient']==patient_temp,])[1] > 0){

    #Append purity to corresponding patient variant record
    snp_whole[i,"Purity_titan"] = unique(bedresult[bedresult['patient']==patient_temp,"titan_purity"])
    #Append local titan call to corresponding patient variant record
    temp_call = bedresult[do.call(paste,bedresult[c("chr_var","start_var","end_var","patient")]) %in% paste(c(paste(snp_whole[i,'Chromosome']), snp_whole[i,'Start_Position'], snp_whole[i,'End_Position'],patient_temp), collapse = " "),"TITAN_call"]
    snp_whole[i,"TITAN_call"] = ifelse(length(temp_call)>0, temp_call, NA)
    #Append local cnv estimates (overall, and strand specific) to corresponding patient variant record
    temp_cnv = bedresult[do.call(paste,bedresult[c("chr_var","start_var","end_var","patient")]) %in% paste(c(paste(snp_whole[i,'Chromosome']), snp_whole[i,'Start_Position'], snp_whole[i,'End_Position'],patient_temp), collapse = " "),c("TITAN_call","local_cnv","minor_cn","major_cn","prevalence")]
    snp_whole[i,"Local.copy.number"] = ifelse(dim(temp_cnv)[1]>0, temp_cnv["local_cnv"], NA)
    snp_whole[i,"minor.cn"] = ifelse(dim(temp_cnv)[1]>0, temp_cnv["minor_cn"], NA) #ifelse(length(temp_cnv)>0, temp_cnv, NA)
    snp_whole[i,"major.cn"] = ifelse(dim(temp_cnv)[1]>0, temp_cnv["major_cn"], NA) # ifelse(length(temp_cnv)>0, temp_cnv, NA)

    #Append local prevalence estimate to corresponding patient variant record
    temp_prev = bedresult[do.call(paste,bedresult[c("chr_var","start_var","end_var","patient")]) %in% paste(c(paste(snp_whole[i,'Chromosome']), snp_whole[i,'Start_Position'], snp_whole[i,'End_Position'],patient_temp), collapse = " "),"prevalence"]
    snp_whole[i,"Prevalence.estimate"] = ifelse(is.na(temp_prev), NA, temp_prev)
    #Append CCF calculation result to corresponding patient variant record (PURITY BASED ON ESTIMATE)
    allele_frac = as.numeric(snp_whole[i,"Variant.Allele.fraction"])
    purity = snp_whole[i,"Purity_titan"]
    if(!(is.na(purity))){if(purity>1){purity=1}}
    copy_num = snp_whole[i,"Local.copy.number"]
    minor_cn = snp_whole[i,"minor.cn"]
    major_cn = snp_whole[i,"major.cn"]
    preval = ifelse(is.na(snp_whole[i,"Prevalence.estimate"]) & (snp_whole[i,"TITAN_call"]=="HET"),1,snp_whole[i,"Prevalence.estimate"])
    snp_whole[i,"Prevalence.estimate"] = preval
    snp_whole[i,"CCF.estimate"] = ifelse(is.na(preval),NA,vaf.to.ccf(allele_frac,purity,preval,minor_cn,major_cn))
    preval=1
    snp_whole[i,"prev1_CCF.estimate"] = ifelse(is.na(preval),NA,vaf.to.ccf(allele_frac,purity,preval,minor_cn,major_cn))

    #Append CCF calculation result to corresponding patient variant record (PURITY BASED ON INPUT)
    allele_frac = as.numeric(snp_whole[i, "Variant.Allele.fraction"])
    purity =  (unique(bedresult[bedresult['patient']==patient_temp, "input_purity"])) #snp_whole[i,"Purity.corrected"]
    if(!(is.na(purity))){if(purity>1){purity=1}}
    snp_whole[i,"input_purity"] = purity
    minor_cn = snp_whole[i,"minor.cn"]
    major_cn = snp_whole[i,"major.cn"]
    preval = ifelse(is.na(snp_whole[i,"Prevalence.estimate"]) & (snp_whole[i,"TITAN_call"]=="HET"),1,snp_whole[i,"Prevalence.estimate"])
    snp_whole[i,"Prevalence.estimate"] = preval
    # cat( " vaf ", allele_frac, " Purity ", purity, " preval ", preval, " minor_cn ", minor_cn, " major_cn ", major_cn, "\n")
    snp_whole[i,"pur.input.CCF.estimate"] = if(is.na(preval) | is.na(purity)){NA}else{vaf.to.ccf(allele_frac,purity,preval,minor_cn,major_cn)}
    preval=1
    snp_whole[i,"pur.input.prev1.CCF.estimate"] = if(is.na(preval) | is.na(purity)){NA}else{vaf.to.ccf(allele_frac,purity,preval,minor_cn,major_cn)}
    
    #Append CCF calculation result to corresponding patient variant record (PURITY BASED ON CORRECTION)
    allele_frac = as.numeric(snp_whole[i, "Variant.Allele.fraction"])
    purity =  (unique(bedresult[bedresult['patient']==patient_temp, "purity_from_mean"])) 
    if(!(is.na(purity))){if(purity>1){purity=1}}
    snp_whole[i,"purity_from_mean"] = purity
    minor_cn = snp_whole[i,"minor.cn"]
    major_cn = snp_whole[i,"major.cn"]
    preval = ifelse(is.na(snp_whole[i,"Prevalence.estimate"]) & (snp_whole[i,"TITAN_call"]=="HET"),1,snp_whole[i,"Prevalence.estimate"])
    snp_whole[i,"Prevalence.estimate"] = preval
    #snp_whole[i,"pur.correc.CCF.estimate"] = ifelse(is.na(preval),NA,vaf.to.ccf(allele_frac,purity,preval,minor_cn,major_cn))
    snp_whole[i,"pur.correc.CCF.estimate"] = if(is.na(preval) | is.na(purity)){NA}else{vaf.to.ccf(allele_frac,purity,preval,minor_cn,major_cn)}
    preval=1
    #snp_whole[i,"pur.correc.prev1.CCF.estimate"] = ifelse(is.na(preval),NA,vaf.to.ccf(allele_frac,purity,preval,minor_cn,major_cn))
    snp_whole[i,"pur.correc.prev1.CCF.estimate"] = if(is.na(preval) | is.na(purity)){NA}else{vaf.to.ccf(allele_frac,purity,preval,minor_cn,major_cn)}
    
    #Append CCF calculation result to corresponding patient variant record BASED ON VAF recalculation AND PURITY BASED ON recalculated VAF
    snp_whole[i,"mean_db_vaf"] = unique(bedresult[bedresult['patient']==patient_temp,"vaf_mean"])
    mean_allele_frac_from_db = unique(bedresult[bedresult['patient']==patient_temp,"vaf_mean"])
    allele_frac = as.numeric(snp_whole[i, "Variant.Allele.fraction"]) #as.numeric(snp_whole[i,"mean_db_vaf"])
    purity = 2* mean_allele_frac_from_db #snp_whole[i,"Purity.estimate"]
    if(!(is.na(purity))){if(purity>1){purity=1}}
    snp_whole[i,"purity.used"] = purity
    snp_whole[i,"vaf.used"] = allele_frac
    minor_cn = snp_whole[i,"minor.cn"]
    major_cn = snp_whole[i,"major.cn"]
    preval = ifelse(is.na(snp_whole[i,"Prevalence.estimate"]) & (snp_whole[i,"TITAN_call"]=="HET"),1,snp_whole[i,"Prevalence.estimate"])
    snp_whole[i,"Prevalence.estimate"] = preval
   # cat( " vaf ", allele_frac, " Purity ", purity, " preval ", preval, " minor_cn ", minor_cn, " major_cn ", major_cn, "\n")
    snp_whole[i,"purity.db.CCF.estimate"] = if(is.na(preval) | is.na(allele_frac) | is.na(purity)){NA}else{vaf.to.ccf(allele_frac,purity,preval,minor_cn,major_cn)}
    preval=1
    snp_whole[i,"purity.db.prev1.CCF.estimate"] = if(is.na(preval) | is.na(allele_frac) | is.na(purity)){NA}else{vaf.to.ccf(allele_frac,purity,preval,minor_cn,major_cn)}
  }
}

write.table(snp_whole,results_output_loh,sep="\t",row.names=F)

#Post preview
##as.data.frame(table(snp[4]))[order(as.data.frame(table(snp[4]))[1]),]
#as.data.frame(table(bedresult['patient']))[order(as.data.frame(table(bedresult['patient']))[1]),]

if(0){
homedir = "/Volumes/Shared/projects/" #Location hosting dlbcl project
home_titan = "dlbcl_relapse_exome/seg_files/TITAN/" #Location in dlbcl project, hosting TITAN output data
localdir = "/Users/jgrewal/Desktop/CCF/"
snp_txt = paste(localdir,"data/snp.txt",sep="")#Variant input file in the format chr start stop patient
snp_complete_txt=paste(localdir,"data/snp_complete.txt",sep="")
bed_output = paste(localdir,"results/local_cn_result.txt",sep="")
results_output = paste(localdir,"results/ccf_result.txt",sep="")
}