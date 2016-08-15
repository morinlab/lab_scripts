#!/usr/bin/env Rscript

# Plot pairwise EXPANDS results for cases with multiple samples
# Takes as input the dm for each sample, performs AF adjustment
# as done in EXPANDS plotting, and them plots pairwise VAF
#
# Mirrors lab_scripts/plot_pyclone_output/plot_ccfs.R
# with EXPANDS-specific modifications

# ------------------ Source functions ---------------------------
command <- commandArgs(trailingOnly = FALSE)
script_name <- sub("--file=", "", command[grep("--file=", command)])
script_path <- dirname(script_name)
functions <- file.path(script_path, "functions.R")
source(functions, chdir = TRUE)


# ------------------- Define arguments --------------------------
p <- arg_parser(description = "Plot EXPANDS results for patients with multiple samples")

# Positional
p <- add_argument(p, arg = "patient", help = "Patient ID")
p <- add_argument(p, arg = "samples", help = "Comma-separated list of sample IDs as used in EXPANDS")
p <- add_argument(p, arg = "dm", help = "Provide dm files for all samples (comma-separated)")
p <- add_argument(p, arg = "mafs",
                  help = "Provide MAFs for all samples (comma-separated, in the same order as samples)")
p <- add_argument(p, arg = "genes", help = "List of genes to label (one per line)")

# Options
p <- add_argument(p, arg = "--out", default = ".",
                  help = "Output directory")
p <- add_argument(p, arg = "--effects", help = "Comma-separated list of VEP effect criteria. If --plot_custom passed and --genes provided, filter mutations to label to the these effects. By default, plots nonsilent variants.",
                  default = "Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,Splice_Site,Translation_Start_Site")

# --------- Get arguments / define other shared variables -------
#args    <- parse_args(p)

args <- parse_args(p, c("FFPE-365", "FFPE-365-A,FFPE-365-R",
                     "6-expands/LOH_and_CN_1_custom/FFPE-365/FFPE-365-A_LOH_state_INDEL_DelMask_maxpm_5_score_2.25_precision_0.05_cnstyle_1.dm.tsv,6-expands/LOH_and_CN_1_custom/FFPE-365/FFPE-365-R_LOH_state_INDEL_DelMask_maxpm_5_score_2.25_precision_0.05_cnstyle_1.dm.tsv",
                     "helper_files/genes.txt"))

patient <- args$patient
samples <- unlist(stri_split(args$samples, regex = ","))
mafs    <- unlist(stri_split(args$mafs, regex = ","))
dm      <- unlist(strsplit(as.character(args$dm), ","))
out_dir <- args$out
effects <- unlist(strsplit(as.character(args$effects), ","))

if (!is.na(args$genes)) {
  genes <- scan(args$genes, what = "character")
} else {
  genes <- NULL
} 

if (length(samples) == 1) {
  print("At least two samples per patient required! Passed:")
  print(samples)
  quit(status = 1)
}

# --------- Plot CCF of all pairs at locus level ----------------
dm_df <- get_loci()

for (i in 1:length(samples)) {
  for (j in 1:length(samples)) {
    
    if (samples[i] == samples[j]) next
    dm_df %>% plot_loci_vaf(c(samples[i], samples[j]))
    
  }
}
