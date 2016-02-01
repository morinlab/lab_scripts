################################################################################
#Setting the working directory
################################################################################
wd <- '/Volumes/morinlab/home/sepideh/exome_coverage/exome_coverage/modifiedcoding'
if(!is.null(wd)) {setwd(wd) }
################################################################################

#Function-1: extract gene coordinates from SureSelect_All_Exon_G3362_with_names_hg19_edited.clean.sorted.bed
################################################################################
fileGene <- 'lymphomapool_coordinates_chr_gene.bed'
fileBed <- 'SureSelect_All_Exon_G3362_with_names_hg19_edited.clean.sorted.bed'
data <- read.table(fileBed, sep='\t', stringsAsFactors=FALSE)
getGeneNames <- function(fileMain, fileGene) {
  names <- read.table(fileGene,  sep="\t")
  names <- levels(names$V4)
  data <- read.table(fileMain, sep='\t', stringsAsFactors=FALSE)
  bed.df <- data.frame(NULL)
  geneTs <-  c(as.character(data$V4))
  listStr <-  strsplit(geneTs, split=':', fixed = TRUE)
  listStr1 <- sapply(listStr, function(x) (x[1]))
  listGene <- strsplit(listStr1, split = 'entg|', fixed = TRUE)
  listGene1 <-  sapply(listGene, function(x) (x[2]))
  rowtokeep <-  is.element(listGene1, names)
  sum(rowtokeep*1)
  out.df <-  data[rowtokeep,]
  out.df[, 4] <-  as.data.frame(listGene1[rowtokeep])
  write.csv(out.df, paste0("lymphoma_",fileMain))
  return(out.df)
}
bed.df <- getGeneNames(fileBed,fileGene)
################################################################################
#a code to find the end and start of gene positions from SureSelect_All_exon...bed 
################################################################################
ag.table <- bed.df
genes <- levels(ag.table$V4)
mins <- c()
maxs <- c()
chrs <- c()
for(gene in genes){
  chr <- unlist(as.character(unique(ag.table[ag.table[,4] == gene, 1])))
  min <- min(ag.table[ag.table[,4] == gene, 2])
  max <- max(ag.table[ag.table[,4] == gene, 3])
  mins <- c(mins,min)
  maxs <- c(maxs,max)
  chrs <- c(chrs,chr)
}
names(mins) <- genes
names(maxs) <- genes
names(chrs) <- genes
bafM.df <- data.frame(matrix(ncol = 3, nrow = length(genes)), row.names = genes)
colnames(bafM.df) <- c('Chr', 'Min', 'Max')
bafM.df$Chr <- chrs
bafM.df$Min <- mins
bafM.df$Max <- maxs
bafM.df
################################################################################
#a code generate a modified lymphoma-coordiane table. SureSelect and Cappseq coordinates vary and this code checks both coordinates to generate the largest gene region
################################################################################
fileGene <- 'lymphomapool_coordinates_chr_gene.bed'
cappseq.df <- read.table(fileGene, sep='\t')
ag.table <- cappseq.df
genes <- levels(ag.table$V4)
mins <- c()
maxs <- c()
chrs <- c()
for(gene in genes){
  chr <- unlist(as.character(unique(ag.table[ag.table[,4] == gene, 1])))
  min <- min(ag.table[ag.table[,4] == gene, 2])
  max <- max(ag.table[ag.table[,4] == gene, 3])
  mins <- c(mins,min)
  maxs <- c(maxs,max)
  chrs <- c(chrs,chr)
}
names(mins) <- genes
names(maxs) <- genes
names(chrs) <- genes
cappseqM.df <- data.frame(matrix(ncol = 3, nrow = length(genes)), row.names = genes)
colnames(cappseqM.df) <- c('Chr', 'Min', 'Max')
cappseqM.df$Chr <- chrs
cappseqM.df$Min <- mins
cappseqM.df$Max <- maxs
cappseqM.df$Chr_2 <- bafM.df$Chr 
cappseqM.df$Min_2 <- bafM.df$Min 
cappseqM.df$Max_2 <- bafM.df$Max 
cappseqM.df
minsM <- c()
maxsM <- c()
chrsM <- c()
cappseqM.df[1, 'Chr_2']
for (rown in rownames(cappseqM.df)){
  chr <-  cappseqM.df[rown, 'Chr_2'] 
  print(chr)
  min <- min(cappseqM.df[rown, c(2,5)])
  print(min)
  max <- max(cappseqM.df[rown, c(3,6)])
  print(max)
  minsM <- c(minsM,min)
  maxsM <- c(maxsM,max)
  chrsM <- c(chrsM,chr)
}
cappexome.df <- data.frame(matrix(ncol = 3, nrow = length(rownames(cappseqM.df))), row.names = rownames(cappseqM.df))
colnames(cappexome.df) <- c('Chr', 'Min', 'Max')
cappexome.df$Chr <- chrsM
cappexome.df$Min <- minsM
cappexome.df$Max <- maxsM
cappexome.df
bafM.df <- cappexome.df
remove <- c("BCL10","EP300","ZFP36L1")
bafM.df <- bafM.df[!row.names(bafM.df)%in%remove,]
bafW.df <- bafM.df
bafM.df 

################################################################################
#codes to parse titan.txt 
################################################################################

#Function-2: calculates avarage of allelic ratios from titan.txt files
geneARAvg <- function(vecAR){
  if (length(vecAR)==0){return(NaN)}
  vecARNoNAN <- vecAR[!is.nan(vecAR)]
  return(mean(vecARNoNAN))
}

#Function-3: extract allelic ratio most distant from 0.5 
geneARMaxDis <- function(vecAR){
  if (length(vecAR)==0){return(NaN)}
  ind <- which.max(abs(vecAR-0.5))
  return(vecAR[ind])
}

#Function-3: calculates avarage of minimum allelic ratios
geneARFlipAvg <- function(vecAR){
  if (length(vecAR)==0){return(NaN)}
  ind <- which(vecAR > 0.5)
  vecAR[ind] <- 1-vecAR[ind]
  vecARNoNAN <- vecAR[!is.nan(vecAR)]
  return(mean(vecARNoNAN))
}

#Function-4: parse titan.txt to calculate different allelic ratio values
patientRatios <- function(patient_file, bafM.df) {
  p.table <- read.csv(patient_file, sep="\t")
  genesARAvg <- c()
  genesARFlipAvg <- c()
  genesARMaxDis <- c()
  titanCalls <- c()
  for (rowbaf in rownames(bafM.df)){
    chrMatch <- (bafM.df[rowbaf, 'Chr'] == p.table$Chr)
    maxMatch <- (bafM.df[rowbaf, 'Max'] >= p.table$Position)
    minMatch <- (bafM.df[rowbaf, 'Min'] <= p.table$Position)
    match <- chrMatch & maxMatch & minMatch
    vecAR <- p.table[match, 'AllelicRatio']
    minBAF <- geneARFlipAvg(vecAR)
    maxDis <- geneARMaxDis(vecAR)
    genesARFlipAvg<- c(genesARFlipAvg, minBAF)
    genesARMaxDis<- c(genesARMaxDis, maxDis)
    genesARAvg<- c(genesARAvg, geneARAvg(vecAR))
    maxDisMatch <- (maxDis == p.table$AllelicRatio) & match
    titan <- as.character(unique(p.table[maxDisMatch, 'TITANcall']))
    if (length(titan) == 0){titan <- NaN}
    titanCalls <- c(titanCalls, titan[1]) #titan[1] is to avoid titalCalls generating more thant two copystate. This can happens when there are maxDisMatch of the same values but different Titancalls.  
  }
  return(list(genesARAvg,genesARFlipAvg,genesARMaxDis,titanCalls))
}
# patient_file = 'PT35T_PT35T_PT35N_PT35N_R0135_cluster_2_titan.txt'
# out <- patientRatios(patient_file, bafM.df)
# genesARAvg <- unlist(out[1])
# genesARFlipAvg <- unlist(out[2])
# genesARMaxDis  <- unlist(out[3])
# titanCalls <- unlist(out[4])
# bafM.df$minBAF.AllelicRatio <- genesARFlipAvg
# bafM.df$TITANcall <- titanCalls
# bafM.df
# bafW.df$Min <- bafM.df$Min - 300000
# bafW.df$Max <- bafM.df$Max + 300000
# outW <- patientRatios(patient_file, bafW.df)
# genesARAvgW <- unlist(outW[1])
# genesARFlipAvgW <- unlist(outW[2])
# genesARMaxDisW <- unlist(outW[3])
# titanCallsW <- unlist(outW[4])
# print(titanCallsW)
# bafW.df$minBAF.AllelicRatioW <- genesARFlipAvgW
# bafW.df$TITANcallW <- titanCallsW
# bafW.df
# nanVec <- is.nan(bafM.df$minBAF.AllelicRatio)
# for (rown in 1:nrow(bafM.df)){
#   if(nanVec[rown]){
#     bafM.df[rown, 'minBAF.AllelicRatioC'] <- bafW.df[rown,'minBAF.AllelicRatioW']
#     bafM.df[rown, 'TITANcallC'] <- bafW.df[rown,'TITANcallW'] 
#   } else {
#     bafM.df[rown, 'minBAF.AllelicRatioC'] <- bafM.df[rown,'minBAF.AllelicRatio']
#     bafM.df[rown, 'TITANcallC'] <- bafM.df[rown,'TITANcall']
#   }
# }
# bafM.df
# bafM.df <- bafM.df[,1:3]
# bafW.df <- bafW.df[,1:3]
bafW.df
################################################################################
#a code to parse titan.txt across a number of patients
################################################################################
files = c('PT003_day15_tumour_na_PT003_normal_na_R0066_cluster_2_titan.txt','PT005T_PT005T_PT005N_PT005N_R0118_cluster_1_titan.txt','PT007T_PT007T_PT007N_PT007N_R0119_cluster_1_titan.txt','PT011T_PT011T_PT011N_PT011N_R0121_cluster_1_titan.txt','PT012T_PT012T_PT012N_PT012N_R0122_cluster_4_titan.txt','PT013T_PT013T_PT013N_PT013N_R0123_cluster_1_titan.txt','PT015_tumour_na_PT015_normal_na_R0069_cluster_3_titan.txt','PT017T_PT017T_PT017N_PT017N_R0124_cluster_3_titan.txt','PT018T_PT018T_PT018N_PT018N_R0125_cluster_2_titan.txt','PT19T_PT19T_PT19N_PT19N_R0126_cluster_3_titan.txt','PT20T_PT20T_PT20N_PT20N_R0127_cluster_2_titan.txt','PT21_tumour_na_PT21_normal_na_R0070_cluster_1_titan.txt','PT22T_PT22T_PT22N_PT22N_R0128_cluster_1_titan.txt','PT23_tumour_na_PT23_normal_na_R0071_cluster_4_titan.txt','PT25T_PT25T_PT25N_PT25N_R0129_cluster_2_titan.txt','PT26T_PT26T_PT26N_PT26N_R0130_cluster_5_titan.txt','PT30T_PT30T_PT30N_PT30N_R0131_cluster_4_titan.txt','PT32T_PT32T_PT32N_PT32N_R0132_cluster_3_titan.txt','PT33T_PT33T_PT33N_PT33N_R0133_cluster_3_titan.txt','PT34T_PT34T_PT34N_PT34N_R0134_cluster_2_titan.txt','PT35T_PT35T_PT35N_PT35N_R0135_cluster_2_titan.txt','PT36T_PT36T_PT36N_PT36N_R0136_cluster_2_titan.txt','PT39T_PT39T_PT39N_PT39N_R0137_cluster_3_titan.txt','PT40T_PT40T_PT40N_PT40N_R0138_cluster_1_titan.txt','PT42_tumour_na_PT42_normal_na_R0072_cluster_1_titan.txt','198T2_198-T2_198-GL_198-GL_R0112_cluster_2_titan.txt','255T1_cluster_3_titan.txt','255T2_cluster_3_titan.txt','260T2_260-T2_260-GL_260-GL_R0113_cluster_2_titan.txt','439T2_439-T2_439-GL_439-GL_R0114_cluster_2_titan.txt')
genenum <- nrow(bafM.df)
n=length(files)
bafPT.df <- data.frame(matrix(ncol= 0, nrow = genenum), row.names = rownames(bafM.df))
for (i in 1:n) {
  out <- patientRatios(files[i], bafM.df)
  genesARAvg <- unlist(out[1])
  genesARFlipAvg <- unlist(out[2])
  genesARMaxDis  <- unlist(out[3])
  titanCalls <- unlist(out[4])
  bafPT.df$minBAF.AllelicRatio <- genesARFlipAvg
  bafPT.df$TITANcall <- titanCalls
  bafPT.df
  #parse titan.txt for expanding gene regions to 300kb. This is to check allelic ratios outside of the gene-regions that lack germ-line SNPs
  bafW.df$Min <- bafM.df$Min - 300000
  bafW.df$Max <- bafM.df$Max + 300000
  outW <- patientRatios(files[i],bafW.df)
  genesARAvgW <- unlist(outW[1])
  genesARFlipAvgW <- unlist(outW[2])
  genesARMaxDisW <- unlist(outW[3])
  titanCallsW <- unlist(outW[4])
  bafPT.df$minBAF.AllelicRatioW <- genesARFlipAvgW
  bafPT.df$TITANcallW <- titanCallsW
  bafPT.df
  nanVec <- is.nan(bafPT.df$minBAF.AllelicRatio)
  for (rown in 1:nrow(bafPT.df)){
    if(nanVec[rown]){
      bafPT.df[rown, 'minBAF.AllelicRatioC'] <- bafPT.df[rown,'minBAF.AllelicRatioW']
      bafPT.df[rown, 'TITANcallC'] <- bafPT.df[rown,'TITANcallW'] 
    } else {
      bafPT.df[rown, 'minBAF.AllelicRatioC'] <- bafPT.df[rown,'minBAF.AllelicRatio']
      bafPT.df[rown, 'TITANcallC'] <- bafPT.df[rown,'TITANcall']
    }
    keeps <- c('minBAF.AllelicRatio', 'TITANcall', 'minBAF.AllelicRatioC', 'TITANcallC')
    bafCc.df <- bafPT.df[,keeps]
    names(bafCc.df) <- c('ARS', 'CNSS', 'ARE', 'CNSE')
  }
  pt <- paste0(sapply(strsplit(files[i], "_"), function(x) (x[1])), '_', names(bafCc.df))
  bafPT.df[,pt] <- bafCc.df
}
bafPT.df
names(bafPT.df)
length(bafPT.df)
keeps <- colnames(bafPT.df[,grepl("_ARS", colnames(bafPT.df))]) 
bafPT_AR.df <- bafPT.df[,keeps]
bafPT_AR.df

keeps <- colnames(bafPT.df[,grepl("_ARE", colnames(bafPT.df))]) 
bafPT_AR_E.df <- bafPT.df[,keeps]
bafPT_AR_E.df

keeps <- colnames(bafPT.df[,grepl("CNSS", colnames(bafPT.df))]) 
bafPT_CNS.df <- bafPT.df[,keeps]
bafPT_CNS.df

keeps <- colnames(bafPT.df[,grepl("CNSE", colnames(bafPT.df))]) 
bafPT_CNS_E.df <- bafPT.df[,keeps]
bafPT_CNS_E.df
################################################################################
# Repeat tumour columns. This is to generate the same color (i.e. copy-number state colors) as the tumour datasets for the plasma samples, for which we miss the titan results.
################################################################################
bafPT_AR_Ex.df <- as.data.frame(rep(bafPT_AR.df[,1:30], times = c(3,2,1,3,2,3,3,3,3,1,3,1,3,1,3,3,1,2,3,1,3,3,5,1,1,1,1,1,1,1)), row.names = row.names(bafPT_AR.df))
length(bafPT_AR_Ex.df)
bafPT_AR_EEX.df <- as.data.frame(rep(bafPT_AR_E.df[,1:30], times = c(3,2,1,3,2,3,3,3,3,1,3,1,3,1,3,3,1,2,3,1,3,3,5,1,1,1,1,1,1,1)), row.names = row.names(bafPT_AR.df)) 
bafPT_CNS_Ex.df  <- as.data.frame(rep(bafPT_CNS.df[,1:30], times = c(3,2,1,3,2,3,3,3,3,1,3,1,3,1,3,3,1,2,3,1,3,3,5,1,1,1,1,1,1,1)), row.names = row.names(bafPT_CNS.df))
bafPT_CNS_EEx.df  <- as.data.frame(rep(bafPT_CNS_E.df[,1:30], times = c(3,2,1,3,2,3,3,3,3,1,3,1,3,1,3,3,1,2,3,1,3,3,5,1,1,1,1,1,1,1)), row.names = row.names(bafPT_CNS.df))
colnames(bafPT_AR_Ex.df) <-   c('03.T15','03.P00','03.P15','05.IT0','05.P15','07.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','15.T00','15.P00','15.P15','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T2','255.T1','255.T2','260.T2','439.T2')
colnames(bafPT_AR_EEX.df) <-  c('03.T15','03.P00','03.P15','05.IT0','05.P15','07.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','15.T00','15.P00','15.P15','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T2','255.T1','255.T2','260.T2','439.T2')
colnames(bafPT_CNS_Ex.df) <-  c('03.T15','03.P00','03.P15','05.IT0','05.P15','07.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','15.T00','15.P00','15.P15','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T2','255.T1','255.T2','260.T2','439.T2')
colnames(bafPT_CNS_EEx.df) <- c('03.T15','03.P00','03.P15','05.IT0','05.P15','07.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','15.T00','15.P00','15.P15','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T2','255.T1','255.T2','260.T2','439.T2')
bafPT_AR_Ex.df
bafPT_AR_EEX.df
bafPT_CNS_Ex.df
bafPT_CNS_EEx.df
################################################################################
## make a dataframe the same size as coverage dataset
################################################################################
codes <- c('02.P00','02.P15','03.T15','03.P00','03.P15','04.IT0','04.P00','04.P15','05.IT0','05.P15','07.T00','09.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','14.P00','14.P15','15.T00','15.P00','15.P15','16.P00','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','27.P00','27.P15','29.P00','29.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','38.P00','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T1','198.T2','255.T1','255.T2','260.T1','260.T2','419.T1','419.T2','439.T1','439.T2','584.T1','584.T2','b01.T0','b02.T0','b02.P','b04.T0','b04.P','b06.T0','b06.P','b12.T0','b12.P','b16.T0','b16.P','b19.T0','b19.P','b22.T0','b22.P')
length(codes)
ptbafPT_AR_Ex.df <- data.frame(matrix(ncol=length(codes), nrow = nrow(bafPT_AR_Ex.df)), row.names = c(rownames(bafPT_AR_Ex.df))) 
colnames(ptbafPT_AR_Ex.df) <- c(codes)
ptbafPT_AR_Ex.df
for (rown in rownames(bafPT_AR_Ex.df)) {
  for (coln in colnames(bafPT_AR_Ex.df)){
    ptbafPT_AR_Ex.df[rown, coln] <- as.character(bafPT_AR_Ex.df[rown, coln])
  }
}
ptbafPT_AR_Ex.df

ptbafPT_AR_EEX.df <- data.frame(matrix(ncol=length(codes), nrow = nrow(bafPT_AR_EEX.df)), row.names = c(rownames(bafPT_AR_EEX.df))) 
colnames(ptbafPT_AR_EEX.df) <- c(codes)
for (rown in rownames(bafPT_AR_EEX.df)) {
  for (coln in colnames(bafPT_AR_EEX.df)){
    ptbafPT_AR_EEX.df[rown, coln] <-as.character(bafPT_AR_EEX.df[rown, coln])
  }
}
ptbafPT_AR_EEX.df

ptbafPT_CNS_Ex.df <- data.frame(matrix(ncol=length(codes), nrow = nrow(bafPT_CNS_Ex.df)), row.names = c(rownames(bafPT_CNS_Ex.df))) 
colnames(ptbafPT_CNS_Ex.df) <- c(codes)
for (rown in rownames(bafPT_CNS_Ex.df)) {
  for (coln in colnames(bafPT_CNS_Ex.df)){
    ptbafPT_CNS_Ex.df[rown, coln] <- as.character(bafPT_CNS_Ex.df[rown, coln])
  }
}
ptbafPT_CNS_Ex.df

ptbafPT_CNS_EEx.df <- data.frame(matrix(ncol=length(codes), nrow = nrow(bafPT_CNS_EEx.df)), row.names = c(rownames(bafPT_CNS_EEx.df))) 
colnames(ptbafPT_CNS_EEx.df) <- c(codes)
for (rown in rownames(bafPT_CNS_EEx.df)) {
  for (coln in colnames(bafPT_CNS_EEx.df)){
    ptbafPT_CNS_EEx.df[rown, coln] <- as.character(bafPT_CNS_EEx.df[rown, coln])
  }
}
ptbafPT_CNS_EEx.df

################################################################################
#a code to parse Titan's seg.txt for the Copy_Number-segment states across differet patients
################################################################################
remove <- c("BCL10","EP300","ZFP36L1")
cappexome.df <- cappexome.df[!row.names(cappexome.df)%in%remove,]
cappexome.df

#Function-5: parse titan.segs to extract titan calls for different genomic segments
patientCNS_short <- function(patient_file, cappexome.df) {
  p.table <- read.csv(patient_file, sep="\t")
  titanCalls_short <- c()
  for (rowbaf in rownames(cappexome.df)){
    print(rowbaf)
    chrMatch <- (cappexome.df[rowbaf, 'Chr'] == p.table$Chromosome)
    maxMatch <- (cappexome.df[rowbaf, 'Max'] <= p.table$End_Position)
    minMatch <- (cappexome.df[rowbaf, 'Min'] >= p.table$Start_Position)
    match <- chrMatch & minMatch & maxMatch
    matchCN <- as.character(unique(p.table[match, 'Copy_Number']))
    matchTC <- as.character(unique(p.table[match, 'TITAN_call'])) 
    if (length(matchCN) == 0){matchCN <- NaN}
    else {
    if (matchCN == '2'& matchTC == 'NLOH'){matchCN <- '50'}
    }
    titanCalls_short <- c(titanCalls_short, matchCN[1])
}
  return(list(titanCalls_short))
}
patient_file = '255T2_na_PT255_N_na_R0099_cluster_5_segs.txt'
p.table <- read.csv(patient_file, sep="\t")
out <- patientCNS_short(patient_file, cappexome.df)
titanCalls_short <- unlist(out)
length(titanCalls_short)
titanCalls_short

################################################################################
#a code to parse Titan's seg.txt for the Copy_Number across differet patients for 300kb larger segments
################################################################################

#Function-6: parse titan.segs to extract CNV calls for different genomic segments expanded 300kb for original-start and -end. This is to maximize the chance of finding gene-regions of interst in the CNV-call segments by titan.  
patientCNS <- function(patient_file, cappexome.df) {
  p.table <- read.csv(patient_file, sep="\t")
  titanCalls <- c()
  for (rowbaf in rownames(cappexome.df)){
    print(rowbaf)
    chrMatch <- (cappexome.df[rowbaf, 'Chr'] == p.table$Chromosome)
    maxMatch <- (cappexome.df[rowbaf, 'Max'] <= p.table$End_Position+300000)
    minMatch <- (cappexome.df[rowbaf, 'Min'] >= p.table$Start_Position-300000)
    match <- chrMatch & minMatch & maxMatch
    matchCN <- as.character(unique(p.table[match, 'Copy_Number']))
    matchTC <- as.character(unique(p.table[match, 'TITAN_call'])) 
    if (length(matchCN) == 0){matchCN <- NaN}
    else {
      if (matchCN == '2'& matchTC == 'NLOH'){matchCN <- '50'}
    }
    titanCalls <- c(titanCalls, matchCN[1])
}
  return(list(titanCalls))
}
patient_file = 'PT35T_PT35T_PT35N_PT35N_R0135_cluster_2_segs.txt'
p.table <- read.csv(patient_file, sep="\t")
out <- patientCNS(patient_file, cappexome.df)
titanCalls <- unlist(out)
titanCalls

files = c('PT003_day15_tumour_na_PT003_normal_na_R0066_cluster_2_segs.txt','PT005T_PT005T_PT005N_PT005N_R0118_cluster_1_segs.txt','PT007T_PT007T_PT007N_PT007N_R0119_cluster_1_segs.txt','PT011T_PT011T_PT011N_PT011N_R0121_cluster_1_segs.txt','PT012T_PT012T_PT012N_PT012N_R0122_cluster_4_segs.txt','PT013T_PT013T_PT013N_PT013N_R0123_cluster_1_segs.txt','PT015_tumour_na_PT015_normal_na_R0069_cluster_3_segs.txt','PT017T_PT017T_PT017N_PT017N_R0124_cluster_3_segs.txt','PT018T_PT018T_PT018N_PT018N_R0125_cluster_2_segs.txt','PT19T_PT19T_PT19N_PT19N_R0126_cluster_3_segs.txt','PT20T_PT20T_PT20N_PT20N_R0127_cluster_2_segs.txt','PT21_tumour_na_PT21_normal_na_R0070_cluster_1_segs.txt','PT22T_PT22T_PT22N_PT22N_R0128_cluster_1_segs.txt','PT23_tumour_na_PT23_normal_na_R0071_cluster_4_segs.txt','PT25T_PT25T_PT25N_PT25N_R0129_cluster_2_segs.txt','PT26T_PT26T_PT26N_PT26N_R0130_cluster_5_segs.txt','PT30T_PT30T_PT30N_PT30N_R0131_cluster_4_segs.txt','PT32T_PT32T_PT32N_PT32N_R0132_cluster_3_segs.txt','PT33T_PT33T_PT33N_PT33N_R0133_cluster_3_segs.txt','PT34T_PT34T_PT34N_PT34N_R0134_cluster_2_segs.txt','PT35T_PT35T_PT35N_PT35N_R0135_cluster_2_segs.txt','PT36T_PT36T_PT36N_PT36N_R0136_cluster_2_segs.txt','PT39T_PT39T_PT39N_PT39N_R0137_cluster_3_segs.txt','PT40T_PT40T_PT40N_PT40N_R0138_cluster_1_segs.txt','PT42_tumour_na_PT42_normal_na_R0072_cluster_1_segs.txt','198T2_198-T2_198-GL_198-GL_R0112_cluster_2_segs.txt','255T1_cluster_3_segs.txt','255T2_cluster_3_segs.txt','260T2_260-T2_260-GL_260-GL_R0113_cluster_2_segs.txt','439T2_439-T2_439-GL_439-GL_R0114_cluster_2_segs.txt')
cns.df <- cappexome.df
genenum <- nrow(cappexome.df)
n <- length(files)
cnsPT.df <- data.frame(matrix(ncol=0, nrow = genenum), row.names = rownames(cns.df))
for (i in 1:n) {
  out <- patientCNS(files[i], cappexome.df)
  titanCalls <- unlist(out)
  cns.df$Copy_Number <- titanCalls
  keeps <- c('Copy_Number')
  cnsM.df <- cns.df[,keeps]
  pt <- paste0(sapply(strsplit(files[i], "_"), function(x) (x[1])), '_CNS')
  cnsPT.df[,pt] <- cnsM.df
}
cnsPT.df
cnsPT.df[cnsPT.df==NaN] <- '100'
cnsPT_Ex.df <- as.data.frame(rep(cnsPT.df[,1:30], times = c(3,2,1,3,2,3,3,3,3,1,3,1,3,1,3,3,1,2,3,1,3,3,5,1,1,1,1,1,1,1)), row.names = row.names(cnsPT.df))
colnames(cnsPT_Ex.df) <- c('03.T15','03.P00','03.P15','05.IT0','05.P15','07.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','15.T00','15.P00','15.P15','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T2','255.T1','255.T2','260.T2','439.T2')
cnsPT_Ex.df
################################################################################
## make a dataframe the same size as cappseqcov dataset
################################################################################
codes <- c('02.P00','02.P15','03.T15','03.P00','03.P15','04.IT0','04.P00','04.P15','05.IT0','05.P15','07.T00','09.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','14.P00','14.P15','15.T00','15.P00','15.P15','16.P00','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','27.P00','27.P15','29.P00','29.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','38.P00','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T1','198.T2','255.T1','255.T2','260.T1','260.T2','419.T1','419.T2','439.T1','439.T2','584.T1','584.T2','b01.T0','b02.T0','b02.P','b04.T0','b04.P','b06.T0','b06.P','b12.T0','b12.P','b16.T0','b16.P','b19.T0','b19.P','b22.T0','b22.P')
length(codes)
cnsPT_EEx.df <- data.frame(matrix(ncol=length(codes), nrow = nrow(cnsPT_Ex.df)), row.names = c(rownames(cnsPT_Ex.df))) 
colnames(cnsPT_EEx.df) <- c(codes)
cnsPT_EEx.df
for (row in rownames(cnsPT_Ex.df)) {
  for (col in colnames(cnsPT_Ex.df)){
    cnsPT_EEx.df[row, col] <- as.character(cnsPT_Ex.df[row, col])
  }
}
ncol(cnsPT_EEx.df)
cnsPT_EEx.df
write.csv(cnsPT_EEx.df, paste0('titan_CNS44','.csv'))
################################################################################
#same titan's copy number state but for actual focal regions not expanded
################################################################################
cnsPT_short.df <- data.frame(matrix(ncol= 0, nrow = genenum), row.names = rownames(cns.df))
for (i in 1:n) {
  out <- patientCNS_short(files[i], cappexome.df)
  titanCalls_short <- unlist(out)
  cns.df$Copy_Number <- titanCalls_short
  keeps <- c('Copy_Number')
  cnsM.df <- cns.df[,keeps]
  pt <- paste0(sapply(strsplit(files[i], "_"), function(x) (x[1])), '_CNS')
  cnsPT_short.df[,pt] <- cnsM.df
}
cnsPT_short.df
cnsPT_short.df[cnsPT_short.df==NaN] <- '100'
cnsPT_Ex_short.df <- as.data.frame(rep(cnsPT.df[,1:30], times = c(3,2,1,3,2,3,3,3,3,1,3,1,3,1,3,3,1,2,3,1,3,3,5,1,1,1,1,1,1,1)), row.names = row.names(cnsPT.df))
colnames(cnsPT_Ex_short.df) <- c('03.T15','03.P00','03.P15','05.IT0','05.P15','07.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','15.T00','15.P00','15.P15','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T2','255.T1','255.T2','260.T2','439.T2')
cnsPT_Ex_short.df
codes <- c('02.P00','02.P15','03.T15','03.P00','03.P15','04.IT0','04.P00','04.P15','05.IT0','05.P15','07.T00','09.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','14.P00','14.P15','15.T00','15.P00','15.P15','16.P00','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','27.P00','27.P15','29.P00','29.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','38.P00','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T1','198.T2','255.T1','255.T2','260.T1','260.T2','419.T1','419.T2','439.T1','439.T2','584.T1','584.T2','b01.T0','b02.T0','b02.P','b04.T0','b04.P','b06.T0','b06.P','b12.T0','b12.P','b16.T0','b16.P','b19.T0','b19.P','b22.T0','b22.P')
length(codes)
cnsPT_EEx_short.df <- data.frame(matrix(ncol=length(codes), nrow = nrow(cnsPT_Ex_short.df)), row.names = c(rownames(cnsPT_Ex_short.df))) 
colnames(cnsPT_EEx_short.df) <- c(codes)
cnsPT_EEx_short.df
for (row in rownames(cnsPT_Ex_short.df)) {
  for (col in colnames(cnsPT_Ex_short.df)){
    cnsPT_EEx_short.df[row, col] <- as.character(cnsPT_Ex_short.df[row, col])
  }
}
ncol(cnsPT_EEx_short.df)
cnsPT_EEx_short.df
length(cnsPT_EEx_short.df)
