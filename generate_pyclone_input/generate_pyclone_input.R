#!/usr/bin/env Rscript

# Generates a PyClone input TSV file.
# Refactors lab_scripts/expands_wrapper/run_expands.R
#
# Usage :
# Rscript generate_pyclone_input.R seg_file seg_input_mode maf_file sample_name output_dir \
#   [--loh {0, 1, 2, 3}] \
#   [--cn_style {1, 2}] \


# ------------------ Source custom EXPANDS utils ----------------
# See http://stackoverflow.com/a/1815743
command <- commandArgs(trailingOnly = FALSE)
script_name <- sub("--file=", "", command[grep("--file=", command)])
script_path <- dirname(script_name)
wrapper_path <- sub("generate_pyclone_input", "expands_wrapper", script_path)
expands_utils <- file.path(wrapper_path, "expandsUtils.R")
source(expands_utils, chdir = TRUE)


# -------------------------- Libraries --------------------------
suppressMessages(suppressWarnings(library(matlab)))
suppressMessages(suppressWarnings(library(argparser)))


# -------------------------- Define arguments -------------------
p <- arg_parser("PyClone Input Generator")

# positional
p <- add_argument(p, "seg", help = "Input segments file")
p <- add_argument(p, "input_mode",
                     help = "Type of seg file: S (Sequenza), I (IGV-friendly seg file), \
                     T (Titan), O (augmented OncoSNP file)")
p <- add_argument(p, "maf", help = "Input MAF file")
p <- add_argument(p, "sample", help = "Sample ID")
p <- add_argument(p, "output_dir",
                     help = "Path to directory where PyClone output will be saved")

# optional
p <- add_argument(p, "--loh", default = 1,
                     help = "0: ignore LOH events, 1: include all copy-neutral LOH segments \
                     and their BAF in clustering; can help resolve clonal clusters with few mutations, (recommended) \
                     2: include deletion LOH only, 3: include all LOH")
p <- add_argument(p, "--cn_style", default = 2,
                  help = "1 for integer values, 2 for rational numbers calculated from CN log ratios (recommended)")

# --------- Get arguments / define other shared variables -------
args <- parse_args(p)

# # for debugging
# args <- parse_args(p, c("../tumour_copy_number/FFPE-121-F_segments.txt", "S",
#                    "../4-clean_maf/FFPE-121-F.clean.maf",
#                    "FFPE", "."))


seg          <- args$seg
input_mode   <- args$input_mode
maf          <- args$maf
sample       <- args$sample
include_loh  <- args$loh
cn_style     <- args$cn_style
out_dir      <- args$output_dir
dir.create(out_dir, recursive = TRUE)

# could be made arguments
min_freq <-  0.1
max_PM <- 5


# ------------------------------------------------------------- 
# Process CNV data from segments file into input CBS matrix for EXPANDS

if (input_mode == "T") {         # Titan
  processed_seg_output <- process_titan_seg(seg, include_loh, cn_style)
} else if (input_mode == "S") {  # Sequenza
  processed_seg_output <- process_sequenza_seg(seg, include_loh, cn_style)
} else if (input_mode == "O") {  # OncoSNP
  processed_seg_output <- process_oncosnp_seg(seg, include_loh, cn_style)
} else if (input_mode == "I") {  # IGV-friendly
  processed_seg_output <- process_igv_seg(seg, include_loh, cn_style)
} else {
  print("Invalid input mode! Choose one of: S, T, I, O. Exiting script.")
  quit(status = 1)
}

seg2 <- processed_seg_output[1]
seg2 <- do.call(rbind, seg2)
loh_snv_data <- processed_seg_output[2]

# Remove? copied and pasted for now
mask_deletions = FALSE

if(cn_style == 2){
  mask_deletions = FALSE
}

if (mask_deletions) {
  seg2.mask = seg2[,4] < 2
  seg2.save = seg2[,"CN_Estimate"]
  
  seg2[seg2.mask,"CN_Estimate"] = 2
  seg2=cbind(seg2,seg2.save)
  colnames(seg2)= c("chr","startpos","endpos","CN_Estimate","CN_Estimate_nomask")
  
}


# -------------------------------------------------------------
# Process SNV data from MAF into input SNV matrix for EXPANDS
# and generate PyClone input file

process_maf_output <- process_maf(maf)
snv_data <- process_maf_output[1]
maf_keep <- process_maf_output[2]

pyclone_input <- generate_pyclone_input(seg, maf_keep, input_mode)

out_pyclone <- paste0(out_dir, "/", sample, "_pyclone_in.tsv")
write.table(pyclone_input, file = out_pyclone,
            sep = "\t", quote = FALSE, row.names = FALSE)

print(paste0("PyClone input written to ", out_pyclone))

