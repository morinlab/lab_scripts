#!/usr/bin/env Rscript
suppressWarnings(suppressPackageStartupMessages(library('argparse')))
suppressWarnings(suppressPackageStartupMessages(library('readr')))
suppressWarnings(suppressPackageStartupMessages(library('dplyr')))
suppressWarnings(suppressPackageStartupMessages(library('tidyr')))
suppressWarnings(suppressPackageStartupMessages(library('magrittr')))
suppressWarnings(suppressPackageStartupMessages(library('ggplot2')))
suppressWarnings(suppressPackageStartupMessages(library('directlabels')))

parser <- ArgumentParser(description='Plot mutations which change in Cancer Cell Fraction (CCF) in consequtive tumour samples.')

parser$add_argument('loci_file', help='Specify PyClone loci file.')
parser$add_argument('maf_file', help='Specify MAF file.')
parser$add_argument('output_file', help='Specify output file for plot.')
parser$add_argument('--plot_width', default = 12, help='Specify output plot width. [%(default)s]')
parser$add_argument('--plot_height', default = 9, help='Specify output plot height. [%(default)s]')
parser$add_argument('--plot_units', default = "in", choices = c("in","cm","mm"),
                    help='Specify output plot units. [%(default)s]')
parser$add_argument('--sample_order', '-s',
                    help='Specified ordering of comma-separated sample names will be used to determine whether a mutations\' CCF change passes fold change threshold. Sample ordering will be preserved in output plots. Defaults to factor sort of "sample_id" column in "loci_file".')
parser$add_argument('--effects', '-c',
                    default = "Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,Splice_Site,Translation_Start_Site",
                    help="Optionally specify MAF variant classifications (comma-separated) to plot. [%(default)s]")
parser$add_argument('--fc_threshold', '-f',
                    type = 'integer',
                    default = 1,
                    help = "Fold change threshold used filter out loci mutations that do not change much between consecutive samples. [%(default)s]")

args <- parser$parse_args()

loci_file <- read_tsv(args$loci_file)
maf_file <- read_tsv(args$maf_file)

if(is.null(args$sample_order)){
  sample_order <- levels(factor(loci_file[,"sample_id"] %>% unlist()))
}else{
  sample_order <- unlist(strsplit(as.character(args$sample_order), ","))
}

effects <- unlist(strsplit(as.character(args$effects), ","))
fc_threshold <- args$fc_threshold
plot_width <- args$plot_width
plot_height <- args$plot_height
plot_units <- args$plot_units
output_file <- args$output_file

# # For testing
# setwd("/Users/prasathpararajalingam/projects/R_scripts/")
# # Rscript plot_ccf_change.R
# #loci_file <- read_tsv("021_loci.tsv")
# #maf_file <- read_tsv("021.maf")
# 
# loci_file <- read_tsv("035_loci.tsv")
# maf_file <- read_tsv("035.maf")
# 
# #sample_order <- c('01-1-021-2_021-BF','01-2-021-1_021-BF','01-3-021-1_021-BF')
# sample_order <- c('01-1-035-2_01-035-BF','01-2-035-2_01-035-BF')
# effects <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation",
#              "Nonsense_Mutation","Nonstop_Mutation","Splice_Site","Translation_Start_Site")
# fc_threshold <- 1
# plot_width <- 12
# plot_height <- 9
# plot_units <- "in"
# output_file <- "./test_plot.png"
# ##################

loci_file <- loci_file %>%
  mutate(sample_id = factor(sample_id, levels = sample_order, ordered=TRUE)) %>%
  group_by(mutation_id) %>%
  arrange(mutation_id, sample_id) %>%
  rename(ccf = cellular_prevalence, ccf_std = cellular_prevalence_std)

maf_file <- maf_file %>%
  select(gene = Hugo_Symbol,
         chr = Chromosome,
         start = Start_Position,
         end = End_Position,
         varclass = Variant_Classification,
         tumour = Tumor_Sample_Barcode,
         normal = Matched_Norm_Sample_Barcode
         ) %>%
  unite(mutation_id, gene, start, sep = "_")

# Identifies loci mutations that display ccf changes above specified
# fold change CCF threshold
mutations_to_keep <- function(loci_mutation, sample_order, fc_threshold){
  for(i in 2:(length(sample_order))){
    lfc <- log2(loci_mutation[i,"ccf"]/loci_mutation[i -1,"ccf"])
    if(lfc >= fc_threshold || lfc <= (-1 * fc_threshold)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
}

mutations <- loci_file %>% select(mutation_id) %>% distinct() %>% unlist()
loci_filtered <- data.frame()
for(mutation in mutations){
  mutations_to_check <- loci_file %>% filter(mutation_id == mutation)
  keep <- mutations_to_keep(mutations_to_check, sample_order, fc_threshold)
  if(keep == TRUE){
    loci_filtered <- loci_filtered %>% bind_rows(mutations_to_check)
  }
}

rm(loci_file)

# Combine filtered loci and filter MAF data frames
loci_maf <- inner_join(maf_file, loci_filtered, by = "mutation_id") %>%
#  separate(mutation_id, c("gene","start"), sep = "_") %>%
  filter(varclass %in% effects) %>%
  select(gene, sample_id, ccf, cluster_id) %>%
  distinct()

# new feature
loci_maf <- inner_join(maf_file, loci_filtered, by = "mutation_id") %>%
  filter(varclass %in% effects)

# count number of genes per cluster and extract genes that have more than one mutation
# in a particular cluster
duplicate_genes <- loci_maf %>%
  group_by(mutation_id, cluster_id) %>%
  summarize(n()) %>%
  ungroup() %>%
  separate(mutation_id, c("gene","start"), sep = "_") %>%
  group_by(cluster_id, gene) %>%
  summarize(gene_count = n()) %>%
  ungroup() %>%
  filter(gene_count > 1) %>%
  select(gene) %>%
  unlist() %>% unname()

# Combine duplicate genes with start site
loci_maf <- loci_maf %>%
  separate(mutation_id, c("gene","start"), sep = "_") %>%
  mutate(gene = ifelse(gene %in% duplicate_genes, paste(gene, start, sep = "_"), gene)) %>%
  select(gene, sample_id, ccf, cluster_id) %>%
  distinct()
  
# Create plot
p <- loci_maf %>%
  ggplot(aes(x = sample_id, y = ccf, group = factor(gene))) +
  geom_line(aes(color = factor(gene))) +
  facet_grid(cluster_id ~ .) +
  theme(legend.position = 'none') +
  geom_dl(aes(label = gene, colour = gene), method = list('last.qp', cex = 0.7, hjust = -0.1)) +
  xlab("Sample ID") +
  ylab("Cancer cell fraction")

# Save plot
ggsave(filename = basename(output_file), plot = p, device = "png",
       path = paste0(dirname(output_file), .Platform$file.sep), 
       width = plot_width, height = plot_height,
       units = plot_units)
