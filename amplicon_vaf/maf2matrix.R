library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

if("-r" %in% args || "--rownames" %in% args){
  variant_list_filename <- ""
  include_rownames_in_output <- TRUE
} else {
  variant_list_filename <- "variant_list.txt"
  include_rownames_in_output <- FALSE
}

if("-c" %in% args || "--coverage" %in% args){
  include_counts_in_output <- TRUE
} else {
  include_counts_in_output <- FALSE
}

weave = function(m1, m2) {
      output_matrix <- matrix(do.call(rbind, list(m1, m2)), nrow = nrow(m1))
      rownames(output_matrix) <- rownames(m1)
      coverage_colnames <- paste0(colnames(m1), rep("_total_coverage", length(colnames(m1))))
      output_colnames <- c(rbind(colnames(m1), coverage_colnames))
      colnames(output_matrix) <- output_colnames
      output_matrix
}

maf <- read.csv(file=args[1], sep="\t", header=T, stringsAsFactors = F)

if(is.null(maf$Variant_ID)) {
  maf$Variant_ID <- paste(as.character(maf$Chromosome), ":", as.character(maf$Start_Position), maf$Reference_Allele, ">", maf$Tumor_Seq_Allele2, sep="")
}

maf$Total_Count <- maf$Reference_Count + maf$Alternate_Count

vaf_matrix <- acast(maf, Variant_ID ~ Sample_ID, value.var = "VAF")

if (include_counts_in_output) {
  total_count_matrix <- acast(maf, Variant_ID ~ Sample_ID, value.var = "Total_Count")
  vaf_matrix <- weave(vaf_matrix, total_count_matrix)
}

if(variant_list_filename != "") {
  write.table(row.names(vaf_matrix), file=variant_list_filename, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
}

write.table(vaf_matrix, quote=F, row.names=include_rownames_in_output, col.names=NA, sep="\t", na="NaN")

