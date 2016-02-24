################################################################################
#Setting my working directory
################################################################################
wd <- '/Volumes/morinlab/home/sepideh/exome_coverage/exome_coverage/modifiedcoding/exome'
if(!is.null(wd)) {setwd(wd)}
files <- list.files(pattern="*.coverage")
################################################################################
#calculate coverage for each gene relative to overall read depth of each sample library (to account for read-coverage varaince due to seq quantity collected for each sample) and log2Ratio Tumour/normal
################################################################################

#Function-1: extract gene coordinates from SureSelect_All_Exon_G3362_with_names_hg19_edited.clean.sorted.bed
################################################################################
getGeneNames <- function(fileMainCov, fileGeneNames) {
  names <- read.table(fileGeneNames,  sep="\t")
  names <- levels(names$V4)
  data <- read.table(fileMainCov, sep='\t', stringsAsFactors=FALSE)
  normal.df <- data.frame(NULL)
  geneTs <-  c(as.character(data$V4))
  listStr <-  strsplit(geneTs, split=':', fixed = TRUE)
  listStr1 <- sapply(listStr, function(x) (x[1]))
  listGene <- strsplit(listStr1, split = 'entg|', fixed = TRUE)
  listGene1 <-  sapply(listGene, function(x) (x[2]))
  rowtokeep <-  is.element(listGene1, names)
  sum(rowtokeep*1)
  out.df <-  data[rowtokeep,]
  out.df[, 4] <-  as.data.frame(listGene1[rowtokeep])
  write.csv(out.df, paste0("ext",fileMainCov))
  return(out.df)
}

#Function-2:Load individual files and calculate the coverage for each gene relative to the mean depth of the library
################################################################################
meanCalc <- function(cov.df) {
  means <- c()
  setGenes <- levels(cov.df[,4])
  for(gene in setGenes){
    sum5 <- sum(cov.df[cov.df[,4] == gene, 5])
    sum6 <- sum(cov.df[cov.df[,4] == gene, 6])
    meanGean <- sum5/sum6
    means <- c(means,meanGean)
  }
  names(means) <- setGenes
  return(means)
}

#Function-3:calculate genome coverage
################################################################################
genomCov <- function(fileMainCov) {
  names <- read.table(fileMainCov, stringsAsFactors=FALSE)
  allMean <-  mean(as.numeric(names$V5))/mean(names$V6)
  return(allMean)
}

################################################################################
## some code that calculates the gene-wise function across a set of samples. The coverage is based on SureSelect.bed file. 
################################################################################
fileGene = 'lymphomapool_coordinates_chr_gene.bed'
names <- read.table(fileGene,  sep="\t")
files_normalCov <- c('PT003N.rmdup.bam.44coverage','PT003N.rmdup.bam.44coverage','PT005N.rmdup.bam.44coverage','PT007N.rmdup.bam.44coverage','PT009N.rmdup.bam.44coverage','PT011N.rmdup.bam.44coverage','PT011N.rmdup.bam.44coverage','PT012N.rmdup.bam.44coverage','PT012N.rmdup.bam.44coverage','PT013N.rmdup.bam.44coverage','PT015N.rmdup.bam.44coverage','PT017N.rmdup.bam.44coverage','PT018N.rmdup.bam.44coverage','PT018N.rmdup.bam.44coverage','PT19N.rmdup.bam.44coverage','PT20N.rmdup.bam.44coverage','PT021N.rmdup.bam.44coverage','PT22N.rmdup.bam.44coverage','PT023N.rmdup.bam.44coverage','PT25N.rmdup.bam.44coverage','PT26N.rmdup.bam.44coverage','PT30N.rmdup.bam.44coverage','PT32N.rmdup.bam.44coverage','PT32N.rmdup.bam.44coverage','PT33N.rmdup.bam.44coverage','PT34N.rmdup.bam.44coverage','PT35N.rmdup.bam.44coverage','PT36N.rmdup.bam.44coverage','PT39N.rmdup.bam.44coverage','PT39N.rmdup.bam.44coverage','PT40N.rmdup.bam.44coverage','PT042N.rmdup.bam.44coverage','198N.srt.rmdup.bam.coverage','255N.srt.rmdup.bam.coverage','255N.srt.rmdup.bam.coverage','260N.srt.rmdup.bam.coverage','419N.srt.rmdup.bam.coverage','419N.srt.rmdup.bam.coverage','439N.srt.rmdup.bam.coverage')
files_tumorCov <-  c('PT003T.rmdup.bam.44coverage','PT003P.rmdup.bam.coverage','PT005T.rmdup.bam.44coverage','PT007T.rmdup.bam.44coverage','PT009T.rmdup.bam.44coverage','PT011T.rmdup.bam.44coverage','PT011P.rmdup.bam.coverage','PT012T.rmdup.bam.44coverage','PT012P.rmdup.bam.coverage','PT013T.rmdup.bam.44coverage','PT015T.rmdup.bam.44coverage','PT017T.rmdup.bam.44coverage','PT018T.rmdup.bam.44coverage','PT18_realign.rmdup.co.bam.coverage','PT19T.rmdup.bam.44coverage','PT20T.rmdup.bam.44coverage','PT021T.rmdup.bam.44coverage','PT22T.rmdup.bam.44coverage','PT023T.rmdup.bam.44coverage','PT25T.rmdup.bam.44coverage','PT26T.rmdup.bam.44coverage','PT30T.rmdup.bam.44coverage','PT32T.rmdup.bam.44coverage','PT032P.rmdup.bam.coverage','PT33T.rmdup.bam.44coverage','PT34T.rmdup.bam.44coverage','PT35T.rmdup.bam.44coverage','PT36T.rmdup.bam.44coverage','PT39T.rmdup.bam.44coverage','PT039P.rmdup.bam.coverage','PT40T.rmdup.bam.44coverage','PT042T.rmdup.bam.44coverage','198T2.srt.rmdup.bam.coverage','255T1.srt.rmdup.bam.coverage','255T2.srt.rmdup.bam.coverage','260T2.srt.rmdup.bam.coverage','419T1.srt.rmdup.bam.coverage','419T2.srt.rmdup.bam.coverage','439T2.srt.rmdup.bam.coverage')	
numcols = length(files_normalCov)
normal.df <- getGeneNames(files_normalCov[38], fileGene) 
normal.df <- normal.df[order(normal.df[,1], normal.df[,2], normal.df[,4]),]# oder functin sorting values first based on column one, then two and then four
meanNormal <- meanCalc(normal.df) # just to count numbers of existing genes
numrows = length(meanNormal)
ptLogRatio.df <- data.frame(matrix(ncol= 0, nrow = numrows), row.names = names(meanNormal))
ptExomCov.df <- data.frame(matrix(ncol= 0, nrow =2 ), row.names =c('normal','tumour'))
for (i in 1:numcols) {
  normal.df <- getGeneNames(files_normalCov[i], fileGene)
  tumor.df <- getGeneNames(files_tumorCov[i], fileGene)
  
  normalGenomCov<- genomCov(files_normalCov[i])
  tumorGenomCov <-  genomCov(files_tumorCov[i])
  
  meanNormal <- meanCalc(normal.df)
  meanTumor <- meanCalc(tumor.df)
  
  corNormal <- meanNormal/normalGenomCov
  print(corNormal)
  corTumor <-  meanTumor/tumorGenomCov
  corRatio <- corTumor/corNormal
  alpha = 1.0
  logRatio <- log2(alpha * corRatio)
  logRatioNoInf <- logRatio[(logRatio != Inf)]
  logRatiomMedian <- logRatio# - median(logRatioNoInf)
  pt <- strsplit(files_tumorCov[i], ".rmdup")[[1]][1]
  ptLogRatio.df[,pt] <- logRatiomMedian
  ptExomCov.df[,pt] <- c(normalGenomCov,tumorGenomCov)
  numrow = length(logRatiomMedian)
}
ptExomCov.df
remove <- c("BCL10","EP300","ZFP36L1") #  Three genes missing for some QCROC samples and therefore are deleted from analysis
ptLogRatio.df <- ptLogRatio.df[!row.names(ptLogRatio.df)%in%remove,]
ptLogRatio.df
write.csv(ptLogRatio.df, paste0('logRatio_tumNor44','.csv'))
################################################################################
## Repeat certain columns
################################################################################
ptLogRatio_Ex.df <- as.data.frame(rep(ptLogRatio.df[,1:39], times = c(1,2,2,1,1,1,2,1,1,3,3,3,1,2,1,3,1,3,1,3,3,1,1,1,3,1,3,3,3,2,1,1,2,1,1,2,1,1,2)), row.names = row.names(ptLogRatio.df))
length(ptLogRatio_Ex.df)
colnames(ptLogRatio_Ex.df) <- c('03.T15','03.P00','03.P15','05.IT0','05.P15','07.T00','09.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','15.T00','15.P00','15.P15','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T1','198.T2','255.T1','255.T2','260.T1','260.T2','419.T1','419.T2','439.T1','439.T2')
ptLogRatio_Ex.df[ptLogRatio_Ex.df==-Inf] <- 'NA'
ptLogRatio_Ex.df
################################################################################
## make a dataframe of exome coverage which includes the same samples and genes as cappseqcov dataset
################################################################################
codes <- c('02.P00','02.P15','03.T15','03.P00','03.P15','04.IT0','04.P00','04.P15','05.IT0','05.P15','07.T00','09.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','14.P00','14.P15','15.T00','15.P00','15.P15','16.P00','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','27.P00','27.P15','29.P00','29.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','38.P00','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T1','198.T2','255.T1','255.T2','260.T1','260.T2','419.T1','419.T2','439.T1','439.T2','584.T1','584.T2','b01.T0','b02.T0','b02.P','b04.T0','b04.P','b06.T0','b06.P','b12.T0','b12.P','b16.T0','b16.P','b19.T0','b19.P','b22.T0','b22.P')
length(codes)
patientcovexome.df <- data.frame(matrix(ncol=length(codes), nrow = nrow(ptLogRatio_Ex.df)), row.names = c(rownames(ptLogRatio_Ex.df))) 
colnames(patientcovexome.df) <- c(codes)
patientcovexome.df
for (rown in rownames(ptLogRatio_Ex.df)) {
  for (coln in colnames(ptLogRatio_Ex.df)){
    patientcovexome.df[rown, coln] <- ptLogRatio_Ex.df[rown, coln]
  }
}
patientcovexome.df
write.csv(patientcovexome.df, paste0('logRatio_tumNor44_Ex','.csv'))
################################################################################
# zscores-coverage _ across patiens for the exome datasets  
################################################################################
data.df <- patientcovexome.df
length(data.df)
colnames(data.df)
rownames(data.df)
length (patientz.df)
colnames(patientz.df)
logRatio_tumNor.df <- data.frame(matrix(ncol=ncol(patientz.df), nrow = nrow(patientz.df)), row.names = c(rownames(patientz.df))) 
colnames(logRatio_tumNor.df) <- c(colnames(patientz.df))
logRatio_tumNor.df 
for (rown in rownames(data.df)) {
  for (coln in colnames(data.df)){
    logRatio_tumNor.df[rown, coln] <- as.character(data.df[rown, coln])
  }
}
logRatio_tumNor.df
################################################################################
# zscoresLogRatio _ across patiens
################################################################################
zLog_tumNor.df <- logRatio_tumNor.df
zLog_tumNor.df
for (i in 1:length(c(rownames(zLog_tumNor.df)))){
  row <- as.numeric(zLog_tumNor.df[i,])
  mean <- mean(row,na.rm=TRUE)
  print(mean)
  standev <- sd(row,na.rm=TRUE)
  print(standev)
  zrow <- (row-mean)/standev
  zrow.df <- t(as.data.frame(zrow)) # when making a datafram for each row it makes a columnwise datafram: zrow with four values framed in four rows that's why we need to transpose the dataframe  
  zLog_tumNor.df[i,] <- zrow.df
}
zLog_tumNor.df
rangecovexome <- range(zLog_tumNor.df[(zLog_tumNor.df != Inf) & !is.na(zLog_tumNor.df)])
rangecovexome
################################################################################
## some code that calculates the exome gene-coverage across a set of samples based on the cappseqcoverage.bed file
################################################################################
files_normal_cappseqCov <- c('PT003N.rmdup.bam.lymphcapcoverage','PT003N.rmdup.bam.lymphcapcoverage','PT005N.rmdup.bam.lymphcapcoverage','PT007N.rmdup.bam.lymphcapcoverage','PT009N.rmdup.bam.lymphcapcoverage','PT011N.rmdup.bam.lymphcapcoverage','PT011N.rmdup.bam.lymphcapcoverage','PT012N.rmdup.bam.lymphcapcoverage','PT012N.rmdup.bam.lymphcapcoverage','PT013N.rmdup.bam.lymphcapcoverage','PT015N.rmdup.bam.lymphcapcoverage','PT017N.rmdup.bam.lymphcapcoverage','PT018N.rmdup.bam.lymphcapcoverage','PT018N.rmdup.bam.lymphcapcoverage','PT19N.rmdup.bam.lymphcapcoverage','PT20N.rmdup.bam.lymphcapcoverage','PT021N.rmdup.bam.lymphcapcoverage','PT22N.rmdup.bam.lymphcapcoverage','PT023N.rmdup.bam.lymphcapcoverage','PT25N.rmdup.bam.lymphcapcoverage','PT26N.rmdup.bam.lymphcapcoverage','PT30N.rmdup.bam.lymphcapcoverage','PT32N.rmdup.bam.lymphcapcoverage','PT32N.rmdup.bam.lymphcapcoverage','PT33N.rmdup.bam.lymphcapcoverage','PT34N.rmdup.bam.lymphcapcoverage','PT35N.rmdup.bam.lymphcapcoverage','PT36N.rmdup.bam.lymphcapcoverage','PT39N.rmdup.bam.lymphcapcoverage','PT39N.rmdup.bam.lymphcapcoverage','PT40N.rmdup.bam.lymphcapcoverage','PT042N.rmdup.bam.lymphcapcoverage','198-N.srt.rmdup.bam.lymphcapcoverage','255-N.srt.rmdup.bam.lymphcapcoverage','255-N.srt.rmdup.bam.lymphcapcoverage','260-N.srt.rmdup.bam.lymphcapcoverage','419-N.srt.rmdup.bam.lymphcapcoverage','419-N.srt.rmdup.bam.lymphcapcoverage','439-N.srt.rmdup.bam.lymphcapcoverage')
files_tumor_cappseqCov <-  c('PT003T.rmdup.bam.lymphcapcoverage','PT003P.rmdup.bam.lymphcapcoverage','PT005T.rmdup.bam.lymphcapcoverage','PT007T.rmdup.bam.lymphcapcoverage','PT009T.rmdup.bam.lymphcapcoverage','PT011T.rmdup.bam.lymphcapcoverage','PT011P.rmdup.bam.lymphcapcoverage','PT012T.rmdup.bam.lymphcapcoverage','PT012P.rmdup.bam.lymphcapcoverage','PT013T.rmdup.bam.lymphcapcoverage','PT015T.rmdup.bam.lymphcapcoverage','PT017T.rmdup.bam.lymphcapcoverage','PT018T.rmdup.bam.lymphcapcoverage','PT018P.rmdup.bam.lymphcapcoverage','PT19T.rmdup.bam.lymphcapcoverage','PT20T.rmdup.bam.lymphcapcoverage','PT021T.rmdup.bam.lymphcapcoverage','PT22T.rmdup.bam.lymphcapcoverage','PT023T.rmdup.bam.lymphcapcoverage','PT25T.rmdup.bam.lymphcapcoverage','PT26T.rmdup.bam.lymphcapcoverage','PT30T.rmdup.bam.lymphcapcoverage','PT32T.rmdup.bam.lymphcapcoverage','PT32P.rmdup.bam.lymphcapcoverage','PT33T.rmdup.bam.lymphcapcoverage','PT34T.rmdup.bam.lymphcapcoverage','PT35T.rmdup.bam.lymphcapcoverage','PT36T.rmdup.bam.lymphcapcoverage','PT39T.rmdup.bam.lymphcapcoverage','PT039P.rmdup.bam.lymphcapcoverage','PT40T.rmdup.bam.lymphcapcoverage','PT042T.rmdup.bam.lymphcapcoverage','198-T2.srt.rmdup.bam.lymphcapcoverage','255-T1.srt.rmdup.bam.lymphcapcoverage','255-T2.srt.rmdup.bam.lymphcapcoverage','260-T2.srt.rmdup.bam.lymphcapcoverage','419-T1.srt.rmdup.bam.lymphcapcoverage','419-T2.srt.rmdup.bam.lymphcapcoverage','439-T2.srt.rmdup.bam.lymphcapcoverage')	
genepool <- get_coverage(files_tumor_cappseqCov[3])
names(genepool)
genenum <- length(genepool)
numcols = length(files_tumor_cappseqCov)
genemat_normal <- matrix(nrow=genenum,ncol=length(files_normal_cappseqCov))
genemat_tumour <- matrix(nrow=genenum,ncol=length(files_tumor_cappseqCov))
rownames(genemat_normal) <- names(genepool)
rownames(genemat_tumour) <- names(genepool)
#gene-coverage across normal exomes
################################################################################
j <- 1
for (i in 1:numcols) {
  normal_cappseq.df <- get_coverage(files_normal_cappseqCov[i])
  namesN <- names(normal_cappseq.df)
    i <- 1
  for (genemeanN in normal_cappseq.df){
    genemat_normal[namesN[i],j] <- genemeanN
    i=i+1
  }
  j=j+1
}
genemat_normal <- as.data.frame(genemat_normal)
genemat_normal
#gene-coverage across tumour exomes
################################################################################
j <- 1
for (i in 1:numcols) {
  tumour_cappseq.df <- get_coverage(files_tumor_cappseqCov[i])
  namesT <- names(tumour_cappseq.df)
    i <- 1
  for (genemeanT in tumour_cappseq.df){
    genemat_tumour[namesT[i],j] <- genemeanT
    i=i+1
  }
  j=j+1
}
genemat_tumour <- as.data.frame(genemat_tumour)
genemat_tumour
#gene log2 coverage-ratio
################################################################################
ptLogRatio_cappcov.df <- data.frame(matrix(ncol=0, nrow = nrow(genemat_tumour)), row.names = c(rownames(genemat_tumour)))
exomeRatio_cappcov <- genemat_tumour/genemat_normal
alpha = 1.0
logRatio_cappcov <- log2(alpha * exomeRatio_cappcov)
logRatioNoInf_cappcov <- logRatio_cappcov[(logRatio_cappcov != Inf)]
logRatiomMedian_cappcov <- logRatio_cappcov# - median(logRatioNoInf_cappseq)
logRatiomMedian_cappcov
for (i in 1:numcols) {
  pt_cappcov <- strsplit(files_tumor_cappseqCov[i], ".rmdup")[[1]][1]
  ptLogRatio_cappcov.df[,pt_cappcov] <- logRatiomMedian_cappcov[,i]
  numrow = length(logRatiomMedian_cappcov)
}
names(ptLogRatio_cappcov.df)
length(ptLogRatio_cappcov.df)
ptLogRatio_cappcov.df
remove <- c("BCL10","EP300","ZFP36L1")
ptLogRatio_cappcov.df <- ptLogRatio_cappcov.df[!row.names(ptLogRatio_cappcov.df)%in%remove,]
write.csv(ptLogRatio_cappcov.df, paste0('logRatio_cappcov4','.csv'))
################################################################################
## Repeat certain columns
################################################################################
ptLogRatio_cappcovEx.df <- as.data.frame(rep(ptLogRatio_cappcov.df[,1:39], times = c(1,2,2,1,1,1,2,1,1,3,3,3,1,2,1,3,1,3,1,3,3,1,1,1,3,1,3,3,3,2,1,1,2,1,1,2,1,1,2)), row.names = row.names(ptLogRatio_cappcov.df))
length(ptLogRatio_cappcovEx.df)
colnames(ptLogRatio_cappcovEx.df) <- c('03.T15','03.P00','03.P15','05.IT0','05.P15','07.T00','09.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','15.T00','15.P00','15.P15','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T1','198.T2','255.T1','255.T2','260.T1','260.T2','419.T1','419.T2','439.T1','439.T2')
ptLogRatio_cappcovEx.df[ptLogRatio_cappcovEx.df==-Inf] <- 'NA'
ptLogRatio_cappcovEx.df
################################################################################
## make a dataframe the same size as cappseqcov dataset
################################################################################
codes <- c('02.P00','02.P15','03.T15','03.P00','03.P15','04.IT0','04.P00','04.P15','05.IT0','05.P15','07.T00','09.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','14.P00','14.P15','15.T00','15.P00','15.P15','16.P00','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','27.P00','27.P15','29.P00','29.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','38.P00','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T1','198.T2','255.T1','255.T2','260.T1','260.T2','419.T1','419.T2','439.T1','439.T2','584.T1','584.T2','b01.T0','b02.T0','b02.P','b04.T0','b04.P','b06.T0','b06.P','b12.T0','b12.P','b16.T0','b16.P','b19.T0','b19.P','b22.T0','b22.P')
length(codes)
patientcovexome_cappcov.df <- data.frame(matrix(ncol=length(codes), nrow = nrow(ptLogRatio_cappcovEx.df)), row.names = c(rownames(ptLogRatio_cappcovEx.df))) 
colnames(patientcovexome_cappcov.df) <- c(codes)
patientcovexome_cappcov.df
for (rown in rownames(ptLogRatio_cappcovEx.df)) {
  for (coln in colnames(ptLogRatio_cappcovEx.df)){
    patientcovexome_cappcov.df[rown, coln] <- ptLogRatio_cappcovEx.df[rown, coln]
  }
}
patientcovexome_cappcov.df
write.csv(patientcovexome_cappcov.df, paste0('logRatio_tumNorcappoc_Ex','.csv'))
################################################################################
# zscores _ across patiens for the exome datasets  
################################################################################
datacappcov.df <- patientcovexome_cappcov.df
logRatiocappcov_tumNor.df <- data.frame(matrix(ncol=ncol(patientcovexome_cappcov.df), nrow = nrow(patientcovexome_cappcov.df)), row.names = c(rownames(patientcovexome_cappcov.df))) 
colnames(logRatiocappcov_tumNor.df) <- c(colnames(patientcovexome_cappcov.df))
length(logRatiocappcov_tumNor.df)
for (rown in rownames(datacappcov.df)) {
  for (coln in colnames(datacappcov.df)){
    logRatiocappcov_tumNor.df[rown, coln] <- as.character(datacappcov.df[rown, coln])
  }
}
logRatiocappcov_tumNor.df
################################################################################
# zscoresLogRatio _ across patiens
################################################################################
zLogcappcov_tumNor.df <- logRatiocappcov_tumNor.df
zLogcappcov_tumNor.df
for (i in 1:length(c(rownames(zLogcappcov_tumNor.df)))){
  row <- as.numeric(zLogcappcov_tumNor.df[i,])
  mean <- mean(row,na.rm=TRUE)
  print(mean)
  standev <- sd(row,na.rm=TRUE)
  print(standev)
  zrow <- (row-mean)/standev
  zrow.df <- t(as.data.frame(zrow))  
  zLogcappcov_tumNor.df[i,] <- zrow.df
}
zLogcappcov_tumNor.df 
rangecovexome <- range(zLogcappcov_tumNor.df[(zLogcappcov_tumNor.df != Inf) & !is.na(zLogcappcov_tumNor.df)])
rangecovexome
