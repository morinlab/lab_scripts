################################################################################
#Setting my working directory
################################################################################
wd <- '/Volumes/morinlab/home/sepideh/exome_coverage/exome_coverage/modifiedcoding'
if(!is.null(wd)) {setwd(wd) }

################################################################################
# extract the lymphoma_genepool coordinates from SureSelect_All_Exon_G3362_with_names_hg19_edited.clean.sorted.bed
################################################################################
fileGene <- 'lymphomapool_coordinates_chr_gene.bed'
fileBed <- 'SureSelect_All_Exon_G3362_with_names_hg19_edited.clean.sorted.bed'
data <- read.table(fileBed, sep='\t', stringsAsFactors=FALSE)

#Function-1:extract the lymphoma_genepool coordinates from SureSelect_All_Exon_G3362_with_names_hg19_edited.clean.sorted.bed
################################################################################
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
#find the end and start of gene positions from SureSelect_All_exon...bed 
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
#Generate a modified lymphoma-gene coordianes. SureSelect and Cappseq coordinates vary and this code checks both coordinates to generate the largest gene region
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
remove <- c("BCL10","EP300","ZFP36L1") # Three genes missing for some QCROC samples and therefore are deleted from analysis
bafM.df <- bafM.df[!row.names(bafM.df)%in%remove,]
bafM.df
################################################################################
# Coverage Analysis
################################################################################

#Function-2:Load individual files and calculate the coverage for each gene relative to the mean depth of the library
################################################################################
get_coverage <- function(bed_file){
  cov.table <- read.table(bed_file, sep="\t")
  genes <- levels(cov.table[,4])
  means <- c()
  allmean <- mean(cov.table[,6])
  for(gene in genes){
    mean <- mean(cov.table[cov.table[,4] == gene, 6])
    meancor <- mean/allmean
    means <- c(means,meancor)
  }
  names(means) <- genes
  return(means)
}

################################################################################
# some code that calculates the gene-wise coverage across a set of samples
################################################################################
files <- list.files(pattern="*.coverage")
files <- c('pt002d00pl_S14.collapsed.bam.coverage','pt02d15pl_S53.collapsed.bam.coverage','pt03d00pl_S01.collapsed.bam.coverage','pt03d15pl_S02.collapsed.bam.coverage','pt004-biopsy-d00_S1.idt44.rmdup.bam.coverage','pt04d00pl_S10.collapsed.bam.coverage','pt04d15pl_S11.collapsed.bam.coverage','pt005-biopsy-d00_S2.idt44.rmdup.bam.coverage','pt05d15pl_S02.collapsed.bam.coverage','pt11d00pl_S15.collapsed.bam.coverage','pt11d15pl_S16.collapsed.bam.coverage','pt012-biopsy-d00_S2.idt44.rmdup.bam.coverage','pt12d00pl_S03.collapsed.bam.coverage','pt13d00pl_S04.collapsed.bam.coverage','pt13d15pl_S05.collapsed.bam.coverage','pt014d00pl_S16.collapsed.bam.coverage','pt14d15pl_S54.collapsed.bam.coverage','pt15d00pl_S06.collapsed.bam.coverage','pt15d15pl_S07.collapsed.bam.coverage','pt016d00pl_S18.collapsed.bam.coverage','pt17d00pl_S13.collapsed.bam.coverage','pt17d15pl_S14.collapsed.bam.coverage','pt018-biopsy-d00_S7.idt44.rmdup.bam.coverage','pt018-plasma-d00_S7.idt44.rmdup.bam.coverage','pt018-plasma-d15_S4.idt44.rmdup.bam.coverage','pt020-biopsy-d00_S4.idt44.rmdup.bam.coverage','pt20d00pl_S15.collapsed.bam.coverage','pt20d15pl_S16.collapsed.bam.coverage','pt021-biopsy-d00_S1.idt44.rmdup.bam.coverage','pt22d00pl_S08.collapsed.bam.coverage','pt22d15pl_S09.collapsed.bam.coverage','pt025-plasma-d00_S10.idt44.rmdup.bam.coverage','pt025-plasma-d15_S11.idt44.rmdup.bam.coverage','pt26d00pl_S10.collapsed.bam.coverage','pt26d15pl_S11.collapsed.bam.coverage','pt027d00pl_S20.collapsed.bam.coverage','pt27d15pl_S55.collapsed.bam.coverage','pt29d00pl_S56.collapsed.bam.coverage','pt029d15pl_S22.collapsed.bam.coverage','pt033-biopsy-d00_S6.idt44.rmdup.bam.coverage','pt33d00pl_S06.collapsed.bam.coverage','pt33d15pl_S07.collapsed.bam.coverage','pt034-biopsy-d00_S5.idt44.rmdup.bam.coverage','pt035-biopsy-d00_S6.idt44.rmdup.bam.coverage','pt35d00pl_S01.collapsed.bam.coverage','pt35d15pl_S02.collapsed.bam.coverage','pt36d00pl_S17.collapsed.bam.coverage','pt36d15pl_S18.collapsed.bam.coverage','pt038d00pl_S24.collapsed.bam.coverage','pt039-biopsy-d00_S6.idt44.rmdup.bam.coverage','pt039-biopsy-d02_S7.idt44.rmdup.bam.coverage','pt039-biopsy-d15_S6.idt44.rmdup.bam.coverage','pt039-plasma-d00_S12.idt44.rmdup.bam.coverage','pt39d15pl_S12.collapsed.bam.coverage','pt198-T1_S13.fixmate.sorted.rmdup.bam.coverage','pt198-T2_S14.fixmate.sorted.rmdup.bam.coverage','pt260-T1_S15.fixmate.sorted.rmdup.bam.coverage','pt260-T2_S16.fixmate.sorted.rmdup.bam.coverage','pt439-T1_S17.fixmate.sorted.rmdup.bam.coverage','pt439-T2_S18.fixmate.sorted.rmdup.bam.coverage','pt584-T1_S19.fixmate.sorted.rmdup.bam.coverage','pt584-T2_S20.fixmate.sorted.rmdup.bam.coverage','bc001_Biopsy_Pre_idt_new_pool.sort.rmdup.bam.coverage','bc002_Biopsy_Pre_idt_new_pool.sort.rmdup.bam.coverage','bc002pl_S7.collapsed.bam.coverage','bc004_Biopsy_Pre_idt_new_pool.sort.rmdup.bam.coverage','bc004pl_S8.collapsed.bam.coverage','bc006_Biopsy_Pre_idt_new_pool.sort.rmdup.bam.coverage','bc006pl_S9.collapsed.bam.coverage','bc012_Biopsy_Pre_idt_new_pool.sort.rmdup.bam.coverage','bc012pl_S2.collapsed.bam.coverage','bc016-tumour-lymphopool-CAPPture_S6.fixmate.sorted.bam.coverage','bc016pl_S3.collapsed.bam.coverage','bc019-tumour-lymphopool-CAPPture_S5.fixmate.sorted.bam.coverage','bc019pl_S2.collapsed.bam.coverage','bc022-tumour-lymphopool-CAPPture_S4.fixmate.sorted.bam.coverage','bc022pl_S1.collapsed.bam.coverage')
genepool <- get_coverage(files[2])
genenum <- length(genepool)
genemat <- matrix(nrow=genenum,ncol=length(files))
rownames(genemat) <- names(genepool)
j <- 1
for (file in files){
  genepool <- get_coverage(file)
  names <- names(genepool)
  i <- 1
  for (genemean in genepool){
    genemat[names[i],j] <- genemean
    i=i+1
  }
  j=j+1
}
genemat.df <- as.data.frame(genemat)
remove <- c("BCL10","EP300","ZFP36L1")
genemat.df <- genemat.df[!row.names(genemat.df)%in%remove,]
nrow(genemat.df)
colnames(genemat.df) <- c('02.P00','02.P15','03.P00','03.P15','04.IT0','04.P00','04.P15','05.IT0','05.P15','11.P00','11.P15','12.IT0','12.P00','13.P00','13.P15','14.P00','14.P15','15.P00','15.P15','16.P00','17.P00','17.P15','18.IT0','18.P00','18.P15','20.IT0','20.P00','20.P15','21.IT0','22.P00','22.P15','25.P00','25.P15','26.P00','26.P15','27.P00','27.P15','29.P00','29.P15','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.P00','36.P15','38.P00','39.IT0','39.IT2','39.I15','39.P00','39.P15','198.T1','198.T2','260.T1','260.T2','439.T1','439.T2','584.T1','584.T2','b01.T0','b02.T0','b02.P','b04.T0','b04.P','b06.T0','b06.P','b12.T0','b12.P','b16.T0','b16.P','b19.T0','b19.P','b22.T0','b22.P')
genemat.df
################################################################################
# zscores _ across patiens for a complete panel of samples. Output NA for samples which there is no cappseq datasets available yet. 
################################################################################
codes <- c('02.P00','02.P15','03.T15','03.P00','03.P15','04.IT0','04.P00','04.P15','05.IT0','05.P15','07.T00','09.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','14.P00','14.P15','15.T00','15.P00','15.P15','16.P00','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','27.P00','27.P15','29.P00','29.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','38.P00','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T1','198.T2','255.T1','255.T2','260.T1','260.T2','419.T1','419.T2','439.T1','439.T2','584.T1','584.T2','b01.T0','b02.T0','b02.P','b04.T0','b04.P','b06.T0','b06.P','b12.T0','b12.P','b16.T0','b16.P','b19.T0','b19.P','b22.T0','b22.P')
length(codes)
patientz.df <- data.frame(matrix(ncol=length(codes), nrow = nrow(genemat.df)), row.names = c(rownames(genemat.df))) 
colnames(patientz.df) <- c(codes)
for (rown in rownames(genemat.df)) {
  for (coln in colnames(genemat.df)){
    patientz.df[rown, coln] <- as.character(genemat.df[rown, coln])
  }
}
nrow(patientz.df)
ncol(patientz.df)
patientcovz.df <- patientz.df
for (i in 1:length(c(rownames(patientcovz.df)))){
  row <- as.numeric(patientcovz.df[i,])
  mean <- mean(row, na.rm=TRUE)
  standev <- sd(row,na.rm=TRUE)
  zrow <- (row-mean)/standev
  zrow.df <- t(as.data.frame(zrow))
  patientcovz.df[i,] <- zrow.df
}
rangecovcappseq <- range(patientcovz.df[(patientcovz.df != Inf) & !is.na(patientcovz.df)])
rangecovcappseq
length(patientcovz.df)
patientcovz.df
