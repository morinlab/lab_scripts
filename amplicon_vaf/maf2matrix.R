library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

maf <- read.csv(file=args[1], sep="\t", header=T, stringsAsFactors = F)

if(is.null(maf$Variant_ID)) {
  maf$Variant_ID <- paste(as.character(maf$Chromosome), ":", as.character(maf$Start_Position), maf$Reference_Allele, ">", maf$Tumor_Seq_Allele2, sep="")
}

vaf_matrix <- acast(maf, Variant_ID ~ Sample_ID, value.var = "VAF")

write.table(vaf_matrix, quote=F, sep="\t")

