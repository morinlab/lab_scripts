#!/usr/bin/env Rscript

# Run EXPANDS and generate a PyClone input TSV file in process
# Refactors lab_scripts/expands_wrapper/run_expands.R
#
# Usage :
# Rscript run_analysis.R seg_file seg_input_mode maf_file sample_name output_dir \
#   [--loh {0, 1, 2, 3}] \
#   [--max_score MAX_SCORE] \
#   [--precision PRECISION] \
#   [--cn_style {1, 2}]

# -------------------------- Libraries --------------------------
library(matlab)
library(expands)
library(argparser)
source("expandsUtils.R")

# -------------------------- Define arguments -------------------
p <- arg_parser("EXPANDS")

# positional
p <- add_argument(p, "seg", help = "Input segments file")
p <- add_argument(p, "input_mode",
                     help = "Type of seg file: S (Sequenza), I (IGV-friendly seg file), \
                     T (Titan), O (augmented OncoSNP file)")
p <- add_argument(p, "maf", help = "Input MAF file")
p <- add_argument(p, "sample", help = "Sample ID")
p <- add_argument(p, "output_dir",
                     help = "Path to directory where all output will be saved")

# optional
p <- add_argument(p, "--loh",
                     help = "0: ignore LOH events, 1: include all copy-neutral LOH segments \
                     and their BAF in clustering; can help resolve clonal clusters with few mutations, \
                     2: include deletion LOH only, 3: include all LOH",
                     default = 1)
p <- add_argument(p, "--max_score", help = "max_score for EXPANDS", default = 2.25)
p <- add_argument(p, "--precision", help = "precision for EXPANDS", default = 0.05)
p <- add_arguemnt(p, "--cn_style", help = "1 for integer values, 2 for rational numbers calculated from CN log ratios")

# --------- Get arguments / define other shared variables -------
# args <- parse_args(p)
args <- parse_args(p, c("~/tmp/DLC_0010.aug.cnvs", "O", "~/tmp/DLC_0010.aug.maf", "DLC_0010",
                        "~/tmp/expands"))
seg         <- args$seg
input_mode  <- args$input_mode
maf         <- args$maf
sample      <- args$sample
include_loh <- args$loh
max_score   <- as.double(args$max_score)
precision   <- as.double(args$precision)
cn_style    <- args$cn_style
out_dir     <- args$output_dir

# could be made arguments
min_freq <-  0.1
max_PM <- 5

# ------------------------------------------------------------- 
# Process CNV data from segments file into input CBS matrix for EXPANDS

list[seg2, loh_snv_data] <- ifelse(input_mode == "T", process_titan_seg(seg, include_loh, cn_style),       # Titan
                            ifelse(input_mode == "S", process_sequenza_seg(seg, include_loh, cn_style),    # Sequenza
                            ifelse(input_mode == "I", process_igv_seg(seg, include_loh, cn_style),         # IGV-friendly
                            ifelse(input_mode == "O", process_oncosnp_seg(seg, include_loh, cn_style)))))  # OncoSNP

# Remove? copied and pasted for now
mask_deletions = 0 

if(cn_style == 2){
  mask_deletions = 0
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
# and genrate PyClone input file

list[snv_data, maf_keep] <- process_maf(maf)
generate_pyclone_input(seg, maf_keep, out_dir, sample)