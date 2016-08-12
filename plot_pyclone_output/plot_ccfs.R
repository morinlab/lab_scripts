#!/usr/bin/env Rscript

# Plot PyClone results for cases with multiple samples

# ------------------ Source functions ---------------------------
command <- commandArgs(trailingOnly = FALSE)
script_name <- sub("--file=", "", command[grep("--file=", command)])
script_path <- dirname(script_name)
functions <- file.path(script_path, "functions.R")
source(functions, chdir = TRUE)


# ------------------- Define arguments --------------------------
p <- arg_parser(description = "Plot PyClone results for patients with multiple samples")

# Positional
p <- add_argument(p, arg = "patient", help = "Patient ID")
p <- add_argument(p, arg = "samples", help = "Comma-separated list of sample IDs as used in PyClone")
p <- add_argument(p, arg = "loci", help = "Path to loci.tsv PyClone output file")
p <- add_argument(p, arg = "genes", help = "List of genes to label (one per line)")

# Options
p <- add_argument(p, arg = "--out", default = ".",
                  help = "Output directory")
p <- add_argument(p, arg = "--min_cluster_size", default = 5,
                  help = "Mirrors PyClone argument: mininum variants required to plot a cluster")
p <- add_argument(p, arg = "--mafs", default = NULL,
                  help = "Provide MAFs for all samples to filter mutations to label by effect (comma-separated, in the same order as samples)")
p <- add_argument(p, arg = "--genes", default = NULL,
                  help = "If --plot_custom passed as flag, label mutations in these genes in custom plot (specify file with one gene per line)")
p <- add_argument(p, arg = "--effects", help = "Comma-separated list of VEP effect criteria. If --plot_custom passed and --genes provided, filter mutations to label to the these effects. By default, plots nonsilent variants.",
                  default = "Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,Splice_Site,Translation_Start_Site")

# --------- Get arguments / define other shared variables -------
args <- parse_args(p)

patient <- args$patient
samples <- unlist(stri_split(args$samples, regex = ","))
loci    <- args$loci 
out_dir <- args$out
min_var <- as.numeric(args$min_cluster_size)
mafs    <- unlist(stri_split(args$mafs, regex = ","))
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
loci_df <- get_loci()
  
for (i in 1:length(samples)) {
  for (j in 1:length(samples)) {
    
    if (samples[i] == samples[j]) next
    plot_loci_ccf(loci_df, c(samples[i], samples[j]))
    
  }
}