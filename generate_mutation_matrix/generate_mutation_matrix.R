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
  -g <genes.txt>, --genes <genes.txt>             List of genes of interest (after considering 
                                                  hotspots and regions).
                                                  Format: Single column, no header.

  -h <hotspots.tsv>, --hotspots <hotspots.tsv>    List of hotspots to consider separately.
                                                  Format: Two columns (gene<tab>codon), no header.

  -c <collapse.txt>, --collapse <collapse.txt>    List of genes for which all mutations should be 
                                                  considered. If --genes is provided, --collapse 
                                                  should represent a subset of --genes. 
                                                  Format: single column, no header.

  -r <regions.bed>, --regions <regions.bed>       List of regions to consider separately. Must use 
                                                  same genome build as MAF file.
                                                  Format: BED file, no header."

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

if (!is.null(args$genes)){
  genes <- fread(args$genes, header = FALSE, col.names = "gene")$gene
}

if (!is.null(args$hotspots)) {
  hotspots <- fread(args$hotspots, col.names = c("Hugo_Symbol", "Codon"),
                    header = FALSE)
}

if (!is.null(args$regions)) {
  regions <- fread(args$regions, header = FALSE, 
                   col.names = c("chrom", "start", "end", "Feature_Name"))
}

if (!is.null(args$collapse)){
  collapse <- fread(args$collapse, header = FALSE, col.names = "gene")$gene
}


# Tidy data ---------------------------------------------------------------

# Restrict to genic mutations (within gene bodies)
maf <- maf[Variant_Classification %in% genic]

# Extract codon information
maf[Variant_Classification == "Missense_Mutation", 
    Codon := sub("p[.][A-Z]+([0-9]+).*", "\\1", HGVSp_Short)]

# Default to not collapsing any gene
maf[, Collapse := FALSE]

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

if (!is.null(args$genes)){
  unannotated_maf <- unannotated_maf[Hugo_Symbol %in% genes]
}

# Annotate genes that should be collapsed
if (!is.null(args$collapse)) {
  unannotated_maf[Hugo_Symbol %in% collapse, Collapse := TRUE]
}


# Annotate remaining mutations --------------------------------------------

# Ignore silent mutations for genes not being collapsed
unannotated_maf <- unannotated_maf[Collapse | Variant_Classification %in% nonsyn]

# Create annotations for collapsed genes
annotated_mafs$collapsed <- 
  unannotated_maf[(Collapse)][, Feature_Name := paste0(Hugo_Symbol, "_All")]
annotated_mafs$collapsed <- annotated_mafs$collapsed[, maf_cols, with = FALSE]

# Create annotations for uncollapsed genes (only nonsynonymous mutations)
annotated_mafs$not_collapsed <- 
  unannotated_maf[(!Collapse)][, Feature_Name := paste0(Hugo_Symbol, "_Nonsyn")]
annotated_mafs$collapsed <- annotated_mafs$collapsed[, maf_cols, with = FALSE]

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
