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
#   [--pyclone_dir PYCLONE_DIR] \
#   [--pyclone_only {TRUE, FALSE}]


# ------------------ Source custom EXPANDS utils ----------------
# See http://stackoverflow.com/a/1815743
command <- commandArgs(trailingOnly = FALSE)
script_name <- sub("--file=", "", command[grep("--file=", command)])
script_path <- dirname(script_name)
expands_utils <- file.path(script_path, "expandsUtils.R")
source(expands_utils, chdir = TRUE)


# -------------------------- Libraries --------------------------
suppressMessages(suppressWarnings(library(matlab)))
suppressMessages(suppressWarnings(library(expands)))
suppressMessages(suppressWarnings(library(argparser)))


# -------------------------- Define arguments -------------------
p <- arg_parser("EXPANDS")

# positional
p <- add_argument(p, "seg", help = "Input segments file")
p <- add_argument(p, "input_mode",
                     help = "Type of seg file: S (Sequenza), I (IGV-friendly seg file), T (Titan), O (augmented OncoSNP file)")
p <- add_argument(p, "maf", help = "Input MAF file")
p <- add_argument(p, "sample", help = "Sample ID")
p <- add_argument(p, "output_dir",
                     help = "Path to directory where all output will be saved")

# optional EXPANDS parameters
p <- add_argument(p, "--loh", default = 1,
                     help = "0: ignore LOH events, 1: include all copy-neutral LOH segments and their BAF in clustering; can help resolve clonal clusters with few mutations, (recommended), 2: include deletion LOH only, 3: include all LOH")
p <- add_argument(p, "--max_score", default = 2.25, help = "max_score for EXPANDS")
p <- add_argument(p, "--precision", default = 0.05, help = "precision for EXPANDS")
p <- add_argument(p, "--cn_style", default = 2,
                  help = "1 for integer values, 2 for rational numbers calculated from CN log ratios (recommended)")

# optional custom plotting parameters
p <- add_argument(p, "--plot_custom", flag = TRUE,
                  help = "Plot cleaner EPXANDS plots (requires dependencies!).")
p <- add_argument(p, "--plot_native", flag = TRUE,
                  help = "Plot native EXPANDS visualizations")
p <- add_argument(p, "--genes", default = NULL,
                  help = "If --plot_custom passed as flag, label mutations in these genes in custom plot (specify file with one gene per line)")
p <- add_argument(p, arg = "--effects", help = "Comma-separated list of VEP effect criteria. If --plot_custom passed and --genes provided, filter mutations to label to the these effects. By default, plots nonsilent variants.",
                  default = "Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,Splice_Site,Translation_Start_Site")
p <- add_argument(p, arg = "--orderBy", help = "If --plot_custom passed as flag, controls ordering of SNVs in custom plot. Options: chr (Chromosome & Start pos), conf (confidence of SP assignment)",
                  default = "conf")


# --------- Get arguments / define other shared variables -------
args <- parse_args(p)

seg          <- args$seg
input_mode   <- args$input_mode
maf          <- args$maf
sample       <- args$sample

include_loh  <- args$loh
max_score    <- as.double(args$max_score)
precision    <- as.double(args$precision)
cn_style     <- args$cn_style
out_dir      <- args$output_dir

plot_native  <- args$plot_native
plot_custom  <- args$plot_custom 
effects      <- unlist(strsplit(as.character(args$effects), ","))
orderBy      <- args$orderBy

if (!is.na(args$genes)) {
  genes <- scan(args$genes, what = "character")
} else {
  genes <- NULL
}

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
print("Completed processing of segments file")

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
print("Completed processing of SNV data")


# -------------------------------------------------------------
# Merge LOH data into SNV matrix

if (include_loh > 0) {
  
  merge_snv <- rbind(loh_snv_data, snv_data)
  loh_string <- "LOH"
  samp_param <- paste0(sample, "_", loh_string, "_state_INDEL_DelMask_maxpm_", max_PM,
                       "_score_", max_score, "_precision_", precision, "_cnstyle_", cn_style)
} else {
  loh_string <- "no_LOH"
  merge_snv <- snv_data
  samp_param <- paste0(sample, "_", loh_string, "_INDEL_DelMask_maxpm_", max_PM,
                      "_score_", max_score, "_precision_", precision, "_cnstyle_", cn_style)
}

merge_snv <- do.call(rbind, merge_snv)


# -------------------------------------------------------------
# Run EXPANDS steps

# 1. Assign each mutation a copy number based on the segmented CN data
# and filter by CN status and AF
dm <- assignQuantityToMutation(merge_snv, seg2, quantityColumnLabel = "CN_Estimate")

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
ii <- which(as.numeric(dm[, "AF_Tumor"]) * as.numeric(dm[, "CN_Estimate"]) < min_freq)
if (length(ii) > 0) {
  print(paste(length(ii), " SNV(s) excluded due to AF*CN below ", min_freq,
              " (SNV can't be explained by an SP present in ",min_freq*100 ,"% or more of the sample)."))
  dm <- dm[-ii, ]
}

if (precision == 1) {
  precision <- 0.1 / log(nrow(dm) / 7) #use default?
  
}

print(paste0("Parameters: ", samp_param))
print("Completed Step 1: Assignment of CN to each mutation")
print("------------------------------------------")


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
print("Completed Step 2: Prediction of sub-populations")
print("------------------------------------------")


# 3. Assign each SNV to one of predicted SPs
aM <- assignMutations(dm, SPs, max_PM = max_PM)

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
print("Completed Step 3: Assignment of SNVs to predicted SPs")
print("------------------------------------------")

# 4. Plot and save results
# Save raw visualization

if (plot_native) {
  out_fig_raw <- paste0(out_dir, "/expands_native_raw_", samp_param, ".pdf")
  pdf(file = out_fig_raw, onefile = FALSE)
  plotSPs(aM$dm, sampleID = sample, cex = 1, rawAF = TRUE)
  dev.off()
  print("Saved raw visualization")

  # Save VAF-corrected visualization
  out_fig <- paste0(out_dir, "/expands_native_adj_", samp_param, ".pdf")
  pdf(file = out_fig, onefile = FALSE)
  plotSPs(aM$dm, sampleID = sample, cex = 1)
  dev.off()
  print("Saved adjusted-AF visualization")

}

if (plot_custom) {
  # Save raw custom visualization
  out_cust_raw <- paste0(out_dir, "/expands_custom_raw_", samp_param, ".pdf")
  pdf(file = out_cust_raw, width = 6, height = 6, onefile = FALSE)
  plot_expands_SPs(aM$dm, sampleID = sample, maf = maf, rawAF = TRUE, orderBy = orderBy, genes = genes, effects = effects)
  dev.off()
  print("Saved raw custom visualization")
  
  # Save VAF-corrected custom visualization
  out_cust <- paste0(out_dir, "/expands_custom_adj_", samp_param, ".pdf")
  dm_string <- paste0(out_dir, "/expands_mutations_adj_", samp_param, ".tsv")
  pdf(file = out_cust, width = 6, height = 6, onefile = FALSE)
  plot_expands_SPs(aM$dm, sampleID = sample, maf = maf, rawAF = FALSE, orderBy = orderBy, genes = genes, effects = effects, dm_string)
  dev.off()
  print("Saved adjusted-AF custom visualization")
  
}

# Save table with mutations assigned to SPs
out_dm <- paste0(out_dir, "/expands_mutations_", samp_param, ".tsv")
suppressWarnings(write.table(aM$dm, file = out_dm, quote = FALSE, sep = "\t", row.names=FALSE))
print(paste0("Mutation details saved to ", out_dm))

# Save details of final SPs
out_final_sps <- paste0(out_dir, "/expands_sps_", samp_param, ".tsv")
write.table(aM$finalSPs, file = out_final_sps, row.names = FALSE, sep = '\t', quote = FALSE)

print(paste0("Final subpopulations predicted in this sample (saved to ", out_final_sps, "):"))
print("Completed Step 4: Visualization of results")
print("------------------------------------------")
