suppressWarnings(library(argparser))
suppressWarnings(library(stringi))
suppressWarnings(library(readr))
suppressWarnings(library(tidyr))
suppressWarnings(library(reshape2))
suppressWarnings(library(ggplot2))
suppressWarnings(library(ggrepel))
suppressWarnings(library(plyr))
suppressWarnings(library(dplyr))


get_loci <- function() {
  
  loci_df <- loci %>%
    read_tsv %>%                                                     # Read PyClone loci.tsv output file
    melt(id.vars = c("mutation_id", "sample_id", "cluster_id")) %>%  # Reshape data into wide format
    dcast(mutation_id + cluster_id ~ sample_id + variable) %>%      
    separate(mutation_id, c("gene", "pos"),
             sep = "_", extra = "warn")                              # Split mutation_id into separate variables
  
  return(loci_df)
}


plot_loci_ccf <- function(loci_df, samples) {
 
  # Columns and labels
  x_ccf <- paste0(samples[1], "_cellular_prevalence")
  xl <- paste0(samples[1], " CCF")
  y_ccf <- paste0(samples[2], "_cellular_prevalence")
  yl <- paste0(samples[2], " CCF")
  loci_df$cluster_id <- factor(loci_df$cluster_id)
  
  # Keep only clusters with n=min_var mutations
  counts <- plyr::count(loci_df, vars = "cluster_id")
  to_keep <- counts[which(counts$freq > min_var), ] %>% 
    dplyr::select(cluster_id) %>% unlist
  loci_df_to_plot <- dplyr::filter(loci_df, cluster_id %in% to_keep)
  
  # Label known genes of interest
  if (is.na(args$mafs)) {  # No MAFs passed, plot all mutations in passed genes
    label_genes <- loci_df %>%
      filter(gene %in% genes)
    
  } else {                 # MAFs provided, restrict genes to label to effects ones
    
    maf_df <- mafs %>% ldply(read_tsv) %>% 
      dplyr::select(Start_Position, Variant_Classification) %>% 
      dplyr::rename(pos = Start_Position) %>% 
      unique()
    
    loci_df$pos <- as.numeric(loci_df$pos)
    label_genes <- loci_df %>% 
      filter(gene %in% genes) %>%
      left_join(maf_df, by = "pos") %>% 
      filter(Variant_Classification %in% effects)
  }
  
  # Plot
  p <- ggplot(loci_df_to_plot, aes_(x = as.name(x_ccf), y = as.name(y_ccf))) +
    geom_point(size = 0.8, aes(colour = cluster_id)) +
    geom_point(data = label_genes,
               aes_(x = as.name(x_ccf), y = as.name(y_ccf), colour = quote(cluster_id))) +
    geom_text_repel(data = label_genes,
                    aes_(x = as.name(x_ccf), y = as.name(y_ccf), label = quote(gene)),
                         nudge_x = 0.05, nudge_y = 0.05) + 
    xlab(xl) + ylab(yl) + labs(colour = "Cluster") +
    coord_fixed(ratio = 1) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1))
  
  out_string <- paste0(samples[2], "vs", samples[1], ".pdf")
  ggsave(p, file = file.path(out_dir, out_string))

}