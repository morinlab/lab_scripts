#!/usr/bin/env Rscript
#
# Usage ---------------------------------------------------------------------------------------
#
# Rscript tx2gene.tsv /path/to/output/prefix /path/to/A01/quant.sf /path/to/A02/quant.sf [...]
#
# The above command will generate two files:
#   - Read counts corrected for size factor: /path/to/output/prefix.norm_counts.tsv
#   - Variance-stabilized expression values: /path/to/output/prefix.vst.tsv
#
# ---------------------------------------------------------------------------------------------


# Load packages -------------------------------------------------------------------------------

# install.packages(c("readr", "tibble"))
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("tximport", "DESeq2"))

suppressPackageStartupMessages({
  library(readr) # Faster than read.delim
  library(tximport)
  library(DESeq2)
  library(tibble)
})


# Parse command-line arguments ----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# Redirect users to usage information if there are no arguments
if (is.null(args[1]) || is.na(args[1])) {
  stop("Run `head process_salmon.R` to see usage information...")
}

paths <- list()
paths$tx2gene <- args[1]
paths$prefix  <- args[2]
paths$ocounts <- paste0(paths$prefix, ".norm_counts.tsv")
paths$ovst    <- paste0(paths$prefix, ".vst.tsv")
paths$quants  <- args[-(1:2)]


# Load input files ----------------------------------------------------------------------------

# Transcript-to-gene map (data frame)
tx2gene <- suppressMessages(read_tsv(file = paths$tx2gene))

# Set names of Salmon quant.sf files to the parent directory name
# Which should be the sample name
names(paths$quants) <- basename(dirname(paths$quants))

# Load in the expression data in the Salmon quant.sf files
message("Loading Salmon quant.sf files...")
txi <- tximport(
  files = paths$quants,
  type = "salmon",
  tx2gene = tx2gene,
  reader = function(...) suppressMessages(read_tsv(..., progress = FALSE)))


# Normalize the expression data ---------------------------------------------------------------

# Create empty data frame with sample metadata
# Necessary for DESeq2
col_data <- data.frame(sample = colnames(txi$counts))

# Load expression count data and estimate size factors
# Using empty colData and design
dds <- DESeqDataSetFromTximport(txi = txi, colData = col_data, design = ~ 1)
dds <- estimateSizeFactors(object = dds)

# Perform variance-stabilizing transformation
tdds <- vst(object = dds)

# Calculate size factor-adjusted counts
counts <- counts(object = dds, normalized = TRUE)


# Output results ------------------------------------------------------------------------------

# Prevent scientific notation
options(scipen = 100)

# Output read counts corrected for size factor
message("Outputting read counts corrected for size factor...")
ocounts <- format(counts, digits = 4)
ocounts <- as.data.frame(ocounts)
ocounts <- rownames_to_column(ocounts, "gene")
write_tsv(x = ocounts, path = paths$ocounts)

# Output variance-stabilized expression values
message("Outputting variance-stabilized expression values...")
ovst <- format(assay(tdds), digits = 4)
ovst <- as.data.frame(ovst)
ovst <- rownames_to_column(ovst, "gene")
write_tsv(x = ovst, path = paths$ovst)
