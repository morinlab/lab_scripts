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
#a code generate a modified lymphoma-coordiane table. SureSelect and Cappseq coordinates vary and I have this code to generates the largest region
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
remove <- c("BCL10","EP300","ZFP36L1")
cappexome.df <- cappexome.df[!row.names(cappexome.df)%in%remove,]
cappexome.df
################################################################################
#a code to parse Sequenza segments.txt for the total Copy_number_tumour across differet patients
################################################################################
genenum <- nrow(cappexome.df)
cnseqPT.df <- data.frame(matrix(ncol= 0, nrow = genenum), row.names = rownames(cappexome.df))

#Function-2: parse segments.txt output from Sequenza to summarize Copy_number_tumour calls of various genes across different patients 
################################################################################
  patientcns <- function(seq_file, cappexome.df) {
  print(seq_file)
  p.table <- read.csv(seq_file, sep="\t")
  seqCalls <- c()
  for (row in rownames(cappexome.df)){
    chrMatch <- (cappexome.df[row, 'Chr'] == p.table$chromosome)
    maxMatch <- (cappexome.df[row, 'Max'] <= p.table$end.pos)
    minMatch <- (cappexome.df[row, 'Min'] >= p.table$start.pos)
    match <- chrMatch & maxMatch & minMatch # This is an array
    matchCNt <- as.character(unique(p.table[match,'CNt']))
    matchA <- na.omit(as.character(unique(p.table[match,'A'])))
    matchB <- na.omit(as.character(unique(p.table[match,'B'])))
    sequenza <- NaN
   if ((length(matchCNt) == 0)| (length(matchA) == 0)| (length(matchB) == 0)){ }
   else {
    if (matchCNt =='2' & matchA == '0' & matchB == '2'){sequenza <- '50'}
    if (matchCNt =='2' & matchA == '2' & matchB == '0'){sequenza <- '50'}
    if (matchCNt =='2' & matchA == '1' & matchB == '1'){sequenza <- '2'}
    if (matchCNt !='2'){sequenza <- matchCNt}
   }
    seqCalls <- c(seqCalls, sequenza)
  }
  return(list(seqCalls))
}

seq_file <- 'PT003_segments.txt'
sequezancns_out <- patientcns(seq_file, cappexome.df)
seqCalls <- unlist(sequezancns_out)
seqCalls

files =c('PT003_segments.txt','PT005_segments.txt', 'PT007_segments.txt', 'PT009_segments.txt','PT011_segments.txt','PT012_segments.txt','PT013_segments.txt', 'PT015_segments.txt', 'PT017_segments.txt', 'PT018_segments.txt', 'PT019_segments.txt', 'PT020_segments.txt', 'PT021_segments.txt', 'PT022_segments.txt', 'PT023_segments.txt', 'PT025_segments.txt', 'PT026_segments.txt', 'PT030_segments.txt', 'PT032_segments.txt', 'PT033_segments.txt', 'PT034_segments.txt', 'PT035_segments.txt', 'PT036_segments.txt', 'PT039_segments.txt', 'PT040_segments.txt', 'PT042_segments.txt','PT198_segments.txt','PT255T1_segments.txt', 'PT255T2_segments.txt', 'PT260_T2_segments.txt', 'PT419_T1_segments.txt', 'PT419_T2_segments.txt','PT439_T2_segments.txt')
length(files)
cns_seq.df <- cappexome.df
cns_seq.df
genenum <- nrow(cappexome.df)
n <- length(files)
cns_seqPT.df <- data.frame(matrix(ncol= 0, nrow = genenum), row.names = rownames(cns_seq.df))
for (i in 1:n) {
  out <- patientcns(files[i], cappexome.df)
  seqCalls <- unlist(out)
  cns_seq.df$CNt <- seqCalls
  keeps <- c('CNt')
  cns_seqM.df <- cns_seq.df[,keeps]
  pt <- paste0(sapply(strsplit(files[i], "segments.txt"), function(x) (x[1])), '_CNt')
  cns_seqPT.df[,pt] <- cns_seqM.df
}
length(cns_seqPT.df)
cns_seqPT.df[cns_seqPT.df==NaN] <- '100'
cns_seqPT.df
cns_seqPT_Ex.df <- as.data.frame(rep(cns_seqPT.df[,1:33], times = c(3,2,1,1,3,2,3,3,3,3,1,3,1,3,1,3,3,1,2,3,1,3,3,5,1,1,2,1,1,2,1,1,2)), row.names = row.names(cns_seqPT.df))
colnames(cns_seqPT_Ex.df) <- c('03.T15','03.P00','03.P15','05.IT0','05.P15','07.T00','09.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','15.T00','15.P00','15.P15','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T1','198.T2','255.T1','255.T2','260.T1','260.T2','419.T1','419.T2','439.T1','439.T2')
cns_seqPT_Ex.df 
codes <- c('02.P00','02.P15','03.T15','03.P00','03.P15','04.IT0','04.P00','04.P15','05.IT0','05.P15','07.T00','09.T00','11.T00','11.P00','11.P15','12.IT0','12.P00','13.T00','13.P00','13.P15','14.P00','14.P15','15.T00','15.P00','15.P15','16.P00','17.T00','17.P00','17.P15','18.IT0','18.P00','18.P15','19.T00','20.IT0','20.P00','20.P15','21.IT0','22.T00','22.P00','22.P15','23.T00','25.T00','25.P00','25.P15','26.T00','26.P00','26.P15','27.P00','27.P15','29.P00','29.P15','30.T00','32.T00','32.P00','33.IT0','33.P00','33.P15','34.IT0','35.IT0','35.P00','35.P15','36.T00','36.P00','36.P15','38.P00','39.IT0','39.IT2','39.I15','39.P00','39.P15','40.T00','42.T00','198.T1','198.T2','255.T1','255.T2','260.T1','260.T2','419.T1','419.T2','439.T1','439.T2','584.T1','584.T2','b01.T0','b02.T0','b02.P','b04.T0','b04.P','b06.T0','b06.P','b12.T0','b12.P','b16.T0','b16.P','b19.T0','b19.P','b22.T0','b22.P')
length(codes)
cns_seqPT_EEx.df <- data.frame(matrix(ncol=length(codes), nrow = nrow(cns_seqPT_Ex.df)), row.names = c(rownames(cns_seqPT_Ex.df))) 
colnames(cns_seqPT_EEx.df) <- c(codes)
cns_seqPT_EEx.df
for (rown in rownames(cns_seqPT_Ex.df)) {
  for (coln in colnames(cns_seqPT_Ex.df)){
    cns_seqPT_EEx.df[rown, coln] <- as.character(cns_seqPT_Ex.df[rown, coln])
  }
}
length(cns_seqPT_EEx.df)
cns_seqPT_EEx.df
write.csv(cns_seqPT_EEx.df, paste0('sequenza_CNt44','.csv'))
