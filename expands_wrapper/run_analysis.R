#!/usr/bin/env Rscript

# Run EXPANDS and generate a PyClone input TSV file in process
# Refactors lab_scripts/expands_wrapper/run_expands.R
#
# Usage :
# Rscript run_analysis.R seg_file seg_input_mode maf_file sample_name output_dir \
#   [--loh {0, 1, 2, 3}] \
#   [--max_score MAX_SCORE] \
#   [--precision PRECISION] \
#   [--cn_style {1, 2}] \
#   [--pyclone_dir PYCLONE_DIR]


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
p <- add_argument(p, "--loh", default = 1,
                     help = "0: ignore LOH events, 1: include all copy-neutral LOH segments \
                     and their BAF in clustering; can help resolve clonal clusters with few mutations, (recommended) \
                     2: include deletion LOH only, 3: include all LOH")
p <- add_argument(p, "--max_score", default = 2.25, help = "max_score for EXPANDS")
p <- add_argument(p, "--precision", default = 0.05, help = "precision for EXPANDS")
p <- add_arguemnt(p, "--cn_style", default = 2,
                  help = "1 for integer values, 2 for rational numbers calculated from CN log ratios (recommended)")
p <- add_argument(p, "--pyclone_dir", default = NULL, help = "Specific output directory for PyClone files")


# --------- Get arguments / define other shared variables -------
args <- parse_args(p)
#args <- parse_args(p, c("~/tmp/DLC_0010.aug.cnvs", "O", "~/tmp/DLC_0010.aug.maf", "DLC_0010",
#                        "~/tmp/expands"))

seg         <- args$seg
input_mode  <- args$input_mode
maf         <- args$maf
sample      <- args$sample
include_loh <- args$loh
max_score   <- as.double(args$max_score)
precision   <- as.double(args$precision)
cn_style    <- args$cn_style
out_dir     <- args$output_dir
pyclone_dir <- ifelse(args$pyclone_dir, args$pyclone_dir, args$output_dir)

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
pyclone_input <- generate_pyclone_input(seg, maf_keep)

out_pyclone <- paste0(out_dir, sample, "_pyclone_in.tsv")
write.table(pyclone_input, file = out_pyclone,
            sep = "\t", quote = FALSE, row.names = FALSE)

print(paste0("PyClone input written to ", out_pyclone))


# -------------------------------------------------------------
# Merge LOH data into SNV matrix

if (include_loh > 0) {
  
  merge_snv <- rbind(loh_snv_data, snv_data)
  #loh_string = "LOH_"
  samp_param <- paste0(sample, loh_string, "_state_INDEL_DelMask_maxpm_", max_PM,
                       "_score_", max_score, ",precision_", precision)
} else {
  loh_string <- "no_LOH"
  merge_snv <- snv_data
  samp_param <- paste0(sample, loh_string, "_noLOH_INDEL_DelMask_maxpm_", max_PM,
                      "_score_", max_score, ",precision_", precision)
}

print("Merged:")
print(merge_snv)
print(".....")
print(seg2)


# -------------------------------------------------------------
# Run EXPANDS steps

# 1. Assign each mutation a copy number based on the segmented CN data
# and filter by CN status and AF
dm <- assignQuantityToMutation(merge_snv, seg2, "CN_Estimate")

plotF <- 1

ii <- which(is.na(dm[, "CN_Estimate"]))
if (length(ii) > 0) {
  print(paste(length(ii), " SNV(s) excluded due to unavailable copy number in that region."))
  dm <- dm[-ii, ]
}
ii <- which(dm[, "CN_Estimate"] < 1)
homDelRegions = c()
if (length(ii) > 0) {
  print(paste(length(ii), " SNV(s) excluded due to homozygous deletions within that region."))
  homDelRegions <- dm[ii, ]
  dm <- dm[-ii, ]
}
ii <- which(dm[, "CN_Estimate"] > max_PM)
if (length(ii) > 0) {
  print(paste(length(ii), " SNV(s) excluded due to high-level amplifications (>", max_PM,
              "copies) within that region. Consider increasing value of parameter max_PM to facilitate inclusion of these SNVs, provided high coverage data (> 150 fold) is available"))
  dm <- dm[-ii, ]
}
ii <- which(dm[, "AF_Tumor"] * dm[, "CN_Estimate"] < min_freq)
if (length(ii) > 0) {
  print(paste(length(ii), " SNV(s) excluded due to AF*CN below ", min_freq,
              " (SNV can't be explained by an SP present in ",min_freq*100 ,"% or more of the sample)."))
  dm <- dm[-ii, ]
}

if (precision == 1) {
  precision <- 0.1 / log(nrow(dm) / 7) #use default?
  
}

print(dm)
print(loh_snv_data)
print(paste0("Parameters: ", samp_param))


# 2. Predict sub-populations
# Compute the cell frequency probability distribution for each mutation
cfd <- computeCellFrequencyDistributions(dm, max_PM, precision, min_CellFreq = min_freq)

# Filter for mutations where cell freq estimation succeeded
toUseIdx <- which(apply(is.finite(cfd$densities), 1, all))

# Cluster cell frequencies
SPs <- clusterCellFrequencies(cfd$densities[toUseIdx, ], precision, min_CellFreq = min_freq)
# Exclude SPs detected at high noise levels
SPs <- SPs[SPs[, "score"] <= max_score, ] 

out_sps <- paste0(out_dir, "/", samp_param, ".unmodified.sps")
write.table(SPs, file = out_sps, row.names = FALSE, sep = '\t', quote = FALSE)

print(paste0("Subpopulations predicted in this sample (saved to ", out_sps, "):"))
print(SPs)


# 3. Assign each SNV to one of predicted SPs
aM <- assignMutations(dm, sps, max_PM = max_PM)

# Unmasks the deletions for vizualization
# can probably be removed
if (mask_deletions) {
  seg2[, "CN_Estimate"] <- seg2[, "CN_Estimate_nomask"]
  dm <- assignQuantityToMutation(merge_snv, seg2, "CN_Estimate")
  some <- dm[, "startpos"] %in% aM$dm[, "startpos"]
  
  aM$dm[, "CN_Estimate"] <- dm[some, "CN_Estimate"]
  
}

# file = "expands_plot_simu_incl_chr6LOH.pdf"
# plot the Expands image showing mutations assigned to their SPs
# aM$dm[,"%maxP"] = 1

# 4. Plot and save results
# Save raw visualization
out_fig_raw <- paste0(out_dir, "/rawPlot_", samp_param, ".pdf")
pdf(out_fig_raw)
plotSPs(aM$dm, sampleID = sample, cex = 1, rawAF = TRUE)
dev.off()

# Save VAF-corrected visualization
out_fig <- paste0(out_dir, "/", samp_param, ".pdf")
pdf(out_fig)
plotSPs(aM$dm, sampleID = sample, cex = 1)
dev.off()

# Save table with mutations assigned to SPs
out_dm <- paste0(out_dir, "/", samp_param, ".dm.tsv")
suppressWarnings(write.table(aM$dm, file = out_dm, quote = FALSE, sep = "\t", row.names=FALSE))
print(paste0("Mutation details saved to ", out_dm))

# Save details of final SPs
out_final_sps <- paste0(out_dir, "/", samp_param, ".final.sps")
write.table(aM$finalSPs, file = out_final_sps, row.names = FALSE, sep = '\t', quote = FALSE)

print(paste0("Final subpopulations predicted in this sample (saved to ", out_final_sps, "):"))
print(aM$finalSPs)