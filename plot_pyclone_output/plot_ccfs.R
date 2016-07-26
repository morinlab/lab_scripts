#!/usr/bin/env Rscript

# Plot PyClone results for cases with multiple samples

# ------------------ Source functions ---------------------------
command <- commandArgs(trailingOnly = FALSE)
script_name <- sub("--file=", "", command[grep("--file=", command)])
script_path <- dirname(script_name)
functions <- file.path(script_path, "plot_ccf_functions.R")
source(functions, chdir = TRUE)


# ------------------- Define arguments --------------------------
p <- arg_parser(description = "Plot PyClone results for patients with multiple samples")

# Positional
p <- add_argument(p, arg = "patient", help = "Patient ID")
p <- add_argument(p, arg = "samples", help = "Comma-separated list of sample IDs as used in PyClone")
p <- add_argument(p, arg = "pyclone_dir", help = "Path to directory containing PyClone results for patient")
p <- add_argument(p, arg = "genes", help = "List of genes to label (one per line)")

# Options
p <- add_argument(p, arg = "--out", default = ".",
                  help = "Output directory")
p <- add_argument(p, arg = "--min_cluster_size", default = 5,
                  help = "Mirrors PyClone argument: mininum variants required to plot a cluster")
p <- add_argument(p, arg = "--mafs", default = NULL,
                  help = "Provide MAFs for all samples to label nonsilent mutations only (comma-separated, in the same order as samples)")


# --------- Get arguments / define other shared variables -------
args <- parse_args(p)

# for debugging
# args <- parse_args(p, c("FFPE-365", "FFPE-365-A,FFPE-365-R",
#                    "../7-pyclone/2-pyclone_working_dir/FFPE-365",
#                    "../helper_files/lymphoma_genes.txt"))
# args$mafs <- "../4-clean_maf/FFPE-365-A.clean.maf,../4-clean_maf/FFPE-365-R.clean.maf"

patient <- args$patient
samples <- unlist(stri_split(args$samples, regex = ","))
pc_dir  <- args$pyclone_dir 
out_dir <- args$out
min_var <- as.numeric(args$min_cluster_size)
genes   <- read_tsv(args$genes) %>% dplyr::select(1) %>% unlist()
mafs    <- unlist(stri_split(args$mafs, regex = ","))

if (length(samples) == 1) {
  print("At least two samples per patient required! Passed:")
  print(samples)
  quit(status = 1)
}

dir.create(out_dir, recursive = TRUE)

nonsilent <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
              "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
              "Splice_Site", "Translation_Start_Site")


# --------- Plot CCF of all pairs at locus level ----------------
loci <- get_loci()
  
for (i in 1:length(samples)) {
  for (j in 1:length(samples)) {
    
    if (samples[i] == samples[j]) next
    plot_loci_ccf(loci, c(samples[i], samples[j]))
    
  }
}
