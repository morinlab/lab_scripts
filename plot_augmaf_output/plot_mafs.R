#!/usr/bin/env Rscript

# Plot VAF vs. VAF for matched samples and optionally,
# plot density, output private gene lists, and output mean/median VAFs 


# ------------------ Source functions ---------------------------
command <- commandArgs(trailingOnly = FALSE)
script_name <- sub("--file=", "", command[grep("--file=", command)])
script_path <- dirname(script_name)
functions <- file.path(script_path, "plot_maf_functions.R")
source(functions, chdir = TRUE)


# ------------------- Define arguments --------------------------
p <- arg_parser(description = "Plot VAFs and other useful results for patients with multiple samples")

# Positional
p <- add_argument(p, arg = "patient", help = "Patient ID")
p <- add_arguemnt(p, arg = "samples", help = "Comma-separated list of sample IDs")
p <- add_argument(p, arg = "mafs", help = "Comma-separated list of augmented MAF paths, in same order as samples")
p <- add_argument(p, arg = "out_dir", help = "Directory to output figures and results")
p <- add_argument(p, arg = "genes", help = "List of genes to label (one per line)")

# Optional
p <- add_argument(p, arg = "--min_vaf", default = 0.05,
                  help = "Variants with VAF > min_vaf in one sample are plot")
p <- add_argument(p, arg = "--min_depth", default = 3, help = "Require at least this depth in normal allele count")
p <- add_argument(p, arg = "--effects", help = "Comma-separated list of VEP effect criteria. Filter mutations to label to the these effects. By default, plots nonsilent variants.",
                  default = "Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,Splice_Site,Translation_Start_Site")

# Flags
p <- add_argument(p, arg = "--density", flag = TRUE, help = "Plot VAF densities")
p <- add_arguemnt(p, arg = "--private", flag = TRUE, help = "Output private gene lists")
p <- add_argument(p, arg = "--vaf_stats", flag = TRUE, help = "Output mean/median VAF for non-private variants in each sample")


# --------- Get arguments / define other shared variables -------
args <- parse_args(p)

patient   <- args$patient
samples   <- args$samples %>% as.character %>% strsplit(",") %>% unlist
mafs      <- args$mafs %>% as.character %>% strsplit(",") %>% unlist
out_dir   <- args$out_dir
genes     <- scan(args$genes, what = "character") %>% unlist
effects   <- args$effects %>% as.character %>% strsplit(",") %>% unlist
min_vaf   <- args$min_vaf
min_depth <- args$min_depth
plot_dens <- args$density
get_priv  <- args$private
vaf_stats <- args$vaf_stats


# --------- Run analysis --------------------------------------

# Plot VAF vs. VAF for all pairs
plot_all_VAF(samples, mafs, min_depth, min_vaf, genes, effects, out_dir)

# Output list of mutations private to each sample
if (get_priv) get_private_mutations(patient, samples, mafs, min_vaf, min_depth, genes) 

# Plot VAF density by sample
if (plot_dens) plot_vaf_density(patient, samples, mafs, min_vaf, min_depth, outdir)

# Output median/mean VAF of non-private mutations for each sample
if (vaf_stats) get_vaf_stats(patient, samples, mafs, min_vaf, min_depth, outdir)











