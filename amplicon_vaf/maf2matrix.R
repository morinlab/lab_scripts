library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

if("-r" %in% args || "--rownames" %in% args){
  variant_list_filename <- ""
  include_rownames_in_output <- TRUE
} else {
  variant_list_filename <- "variant_list.txt"
  include_rownames_in_output <- FALSE
}

maf <- read.csv(file=args[1], sep="\t", header=T, stringsAsFactors = F)

if(is.null(maf$Variant_ID)) {
  maf$Variant_ID <- paste(as.character(maf$Chromosome), ":", as.character(maf$Start_Position), maf$Reference_Allele, ">", maf$Tumor_Seq_Allele2, sep="")
}

vaf_matrix <- acast(maf, Variant_ID ~ Sample_ID, value.var = "VAF")

if(variant_list_filename != "") {
  write.table(row.names(vaf_matrix), file=variant_list_filename, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
}

write.table(vaf_matrix, quote=F, row.names=include_rownames_in_output, sep="\t", na="NaN")

