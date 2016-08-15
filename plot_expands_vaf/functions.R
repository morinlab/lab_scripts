suppressWarnings(library(argparser))
suppressWarnings(library(stringi))
suppressWarnings(library(readr))
suppressWarnings(library(tidyr))
suppressWarnings(library(reshape2))
suppressWarnings(library(ggplot2))
suppressWarnings(library(ggrepel))
suppressWarnings(library(plyr))
suppressWarnings(library(dplyr))


get_dm <- function() {
  
  # Read all MAFs
  maf_df <- mafs %>% ldply(read_tsv) %>% 
    dplyr::select(Tumor_Sample_Barcode, Start_Position, Hugo_Symbol)
    
  # Read in all dm files into one data frame
  dm_df <- dm %>%
    {
      names(.) <- samples[seq_along(samples)]
      return(.)
    } %>% 
    ldply(read_tsv, .id = "Sample") %>%
    dplyr::select(chr, startpos, Sample, AF_Tumor_Adjusted) %>% 
    mutate(Sample = paste0(Sample, ".VAF")) %>% 
    spread(Sample, AF_Tumor_Adjusted) %>%
    left_join(maf_df, by = c("startpos" = "Start_Position")) %>%
    filter(!is.na(Hugo_Symbol))
  
  return(dm_df)
}

plot_loci_vaf <- function(dm_df, samples) {
 
  # Columns and labels
  x_vaf <- paste0(samples[1], ".VAF")
  xl <- paste0(samples[1], " Adjusted VAF")
  y_vaf <- paste0(samples[2], ".VAF")
  yl <- paste0(samples[2], " Adjusted VAF")
  
  # Label known genes of interest
  maf_df <- mafs %>% ldply(read_tsv) %>% 
      dplyr::select(Start_Position, Variant_Classification) %>% 
      unique()
  
  dm_df$startpos <- as.numeric(dm_df$startpos)
  label_genes <- dm_df %>% 
    filter(Hugo_Symbol %in% genes) %>%
    left_join(maf_df, by = c("startpos" = "Start_Position")) %>% 
    filter(Variant_Classification %in% effects)
  
  # Plot
  p <- ggplot(dm_df, aes_(x = as.name(x_vaf), y = as.name(y_vaf))) +
    geom_point(size = 0.8) +
    geom_point(data = label_genes,
               aes_(x = as.name(x_vaf), y = as.name(y_vaf))) +
    geom_text_repel(data = label_genes,
                    aes_(x = as.name(x_vaf), y = as.name(y_vaf), label = quote(Hugo_Symbol)),
                         nudge_x = 0.05, nudge_y = 0.05) + 
    xlab(xl) + ylab(yl) + coord_fixed(ratio = 1) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1))
  
  out_string <- paste0(samples[2], "vs", samples[1], ".pdf")
  ggsave(p, file = file.path(out_dir, out_string))

}