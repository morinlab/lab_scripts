#!/usr/bin/env Rscript
# Author: Bruno Grande (bgrande@sfu.ca)


# Load libraries ----------------------------------------------------------

suppressWarnings(suppressPackageStartupMessages({
  library(docopt)
  library(data.table)
}))


# Parse arguments ---------------------------------------------------------

ui <- "
Usage: 
  generate_mutation_matrix.R [options] <maf> <output>

Options:
  -h <hotspots.tsv>, --hotspots <hotspots.tsv>    Mutation hotspots to output as separate rows in
                                                  the generated mutation matrix.
                                                  Format: Two columns (gene<tab>codon), no header.
                                                  Default: No hotspots.
  
  -r <regions.bed>, --regions <regions.bed>       Regions (using same genome build) to output as 
                                                  separate rows in the generated mutation matrix.
                                                  Format: BED file with name column, no header.
                                                  Default: No regions.

  -n <nonsyn.txt>, --nonsyn <nonsyn.txt>          Genes of interest for which only nonsynonymous 
                                                  mutations are considered (after removing any
                                                  mutations in the specified hotspots or regions).
                                                  Format: Single column, no header.
                                                  Default: All genes listed in the input MAF file.

  -a <all.txt>, --all <all.txt>                   Genes of interest for which all mutations are 
                                                  considered (after removing any mutations in the 
                                                  specified hotspots or regions), including
                                                  synonymous, nonsynonymous and flanking variants.
                                                  This option takes precedence for genes listed for
                                                  both the --nonsyn and --all options.
                                                  Format: Single column, no header.
                                                  Default: No genes."

args <- docopt(ui)


# Constants ---------------------------------------------------------------

genic <- 
  c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", 
    "Missense_Mutation", "Nonsense_Mutation", "Silent", "Splice_Site", 
    "Translation_Start_Site", "Nonstop_Mutation", "3'UTR", "5'UTR", 
    "3'Flank", "5'Flank", "Intron")

silent <- c("3'UTR", "5'UTR", "3'Flank", "5'Flank", "Intron", "Silent")

nonsyn <- setdiff(genic, silent)


# Load data ---------------------------------------------------------------

maf <- fread(args$maf)

if (!is.null(args$nonsyn)){
  nonsyn <- fread(args$nonsyn, header = FALSE, col.names = "gene")$gene
}

if (!is.null(args$hotspots)) {
  hotspots <- fread(args$hotspots, col.names = c("Hugo_Symbol", "Codon"),
                    header = FALSE)
}

if (!is.null(args$regions)) {
  regions <- fread(args$regions, header = FALSE, 
                   col.names = c("chrom", "start", "end", "Feature_Name"))
}

if (!is.null(args$all)){
  all <- fread(args$all, header = FALSE, col.names = "gene")$gene
} else {
  all <- vector("character")
}


# Tidy data ---------------------------------------------------------------

# Restrict to genic mutations (within gene bodies)
maf <- maf[Variant_Classification %in% genic]

# Extract codon information
maf[Variant_Classification == "Missense_Mutation", 
    Codon := sub("p[.][A-Z]+([0-9]+).*", "\\1", HGVSp_Short)]

# Default to not collapsing any gene
maf[, All := FALSE]

# Set unique key for MAF file
maf[, Mutation_ID := .I]
setkey(maf, Mutation_ID)


# Annotate special mutations ----------------------------------------------

annotated_mafs <- list()

unannotated_maf <- maf

maf_cols <- c(colnames(maf), "Feature_Name")

# Annotate hotspots
if (!is.null(args$hotspots)) {
  hotspots[, Codon := as.character(Codon)]
  hotspots[, Feature_Name := paste0(Hugo_Symbol, "_Codon", Codon)]
  annotated_mafs$hotspots <- merge(unannotated_maf, hotspots, 
                                   by = c("Hugo_Symbol", "Codon"))
  annotated_mafs$hotspots <- annotated_mafs$hotspots[, maf_cols, with = FALSE]
  # Prevent hotspot mutations from being considered elsewhere
  setkey(annotated_mafs$hotspots, Mutation_ID)
  unannotated_maf <- unannotated_maf[!annotated_mafs$hotspots]
}

# Annotate custom regions
if (!is.null(args$regions)) {
  # Switch from 0-based to 1-based indexing
  regions[, start := start + 1]
  regions[, chrom := as.character(chrom)]
  setkey(regions, chrom, start, end)
  annotated_mafs$regions <- foverlaps(
    unannotated_maf, regions, nomatch = 0, 
    by.x = c("Chromosome", "Start_Position", "End_Position"))
  annotated_mafs$regions <- annotated_mafs$regions[, maf_cols, with = FALSE]
  # Prevent hotspot mutations from being considered elsewhere
  setkey(annotated_mafs$regions, Mutation_ID)
  unannotated_maf <- unannotated_maf[!annotated_mafs$regions]
  # Ensure that mutations aren't being picked up more than once
  if (annotated_mafs$regions[, any(duplicated(Mutation_ID))]) {
    stop(paste0(
      "Some mutations were captured by more than one region. Please ensure that ",
      "your regions are disjoint (e.g., by merging them beforehand)."))
  }
}

if (!is.null(args$nonsyn)){
  unannotated_maf <- unannotated_maf[Hugo_Symbol %in% c(nonsyn, all)]
}

# Annotate genes that should be collapsed
if (!is.null(args$all)) {
  unannotated_maf[Hugo_Symbol %in% all, All := TRUE]
}


# Annotate remaining mutations --------------------------------------------

# Ignore silent mutations for genes not being collapsed
unannotated_maf <- unannotated_maf[All | Variant_Classification %in% nonsyn]

# Create annotations for collapsed genes
annotated_mafs$all <- 
  unannotated_maf[(All)][, Feature_Name := paste0(Hugo_Symbol, "_All")]
annotated_mafs$all <- annotated_mafs$all[, maf_cols, with = FALSE]

# Create annotations for uncollapsed genes (only nonsynonymous mutations)
annotated_mafs$nonsyn <- 
  unannotated_maf[(!All)][, Feature_Name := paste0(Hugo_Symbol, "_Nonsyn")]
annotated_mafs$nonsyn <- annotated_mafs$nonsyn[, maf_cols, with = FALSE]

# Merge everything once again (including collapsed and non-collapsed genes)
annotated_maf <- rbindlist(annotated_mafs)

# Check that each mutation only appears once
if (any(annotated_maf[, duplicated(Mutation_ID)])) {
  stop("This error shouldn't be happening. Please report this to the developer.")
}


# Output mutation matrix --------------------------------------------------

# Create matrix
mutation_matrix <- 
  dcast.data.table(annotated_maf, Feature_Name ~ Tumor_Sample_Barcode,
                   fun.aggregate = length, value.var = "Mutation_ID")

fwrite(mutation_matrix, args$output, sep = "\t")
