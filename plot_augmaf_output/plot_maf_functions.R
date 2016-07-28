suppressWarnings(library(argparser))
suppressWarnings(library(stringi))
suppressWarnings(library(readr))
suppressWarnings(library(ggplot2))
suppressWarnings(library(ggrepel))

read_maf <- function(maf_path, min_depth) {
  # Read in a MAF from path, calculate VAF, and threshold on min_depth
  
  maf_df <- read_tsv(maf_path) %>% 
    dplyr::mutate(VAF = ifelse(t_depth > 0, t_alt_count/t_depth, 0)) %>% 
    dplyr::filter(n_depth > min_depth)
  
  return(maf_df)
}

get_maf_overlap <- function(maf_1, maf_2, min_vaf) {
  # Get the overlap between two MAFs given maf_df's
  
  # Get overlap and threshold
  overlap <- dplyr::full_join(maf_1, maf_2,
                              by = c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
                                     "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification")) %>% 
             dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification,
                           t_ref_count.x, t_alt_count.x, VAF.x, t_ref_count.y, t_alt_count.y, VAF.y)
  overlap$maxVAF <- apply(overlap[, c("VAF.x", "VAF.y")], 1, max)
  overlap %<>% dplyr::filter(maxVAF >= min_vaf) 
  
  return(overlap)
}

plot_VAF_vs_VAF <- function(variants, patient, id_1, id_2, genes, effects, out_dir) {
  # Generate and save VAF vs. VAF plots given a variants df
  # returned from get_maf_overlap and the sample IDs
  
  # Scatterplot
  xl = paste0(id_1, " VAF")
  yl = paste0(id_2, " VAF")
  variants$VAF.x <- as.numeric(variants$VAF.x)
  variants$VAF.y <- as.numeric(variants$VAF.y)
  
  label_genes <- variants %>%
    dplyr::filter(Variant_Classification %in% effects) %>% 
    dplyr::filter(Hugo_Symbol %in% genes)
  
  p <- ggplot(variants, aes(x = VAF.x, y = VAF.y)) +
    geom_point(size = 0.8, colour = "gray63") +
    geom_point(data = label_genes, aes(x = VAF.x, y = VAF.y), colour = "red") +
    geom_text_repel(data = label_genes,
                    aes(x = VAF.x, y = VAF.y, label = Hugo_Symbol), color = "blue") +
    labs(title = patient) + xlab(xl) + ylab(yl) + coord_fixed(ratio = 1) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1))
  
  ggsave(p, file = paste0(out_dir, "/", id_2, "vs", id_1, ".png"))
  
  # Smooth scatterplot
  png(file = paste0(out_dir, "/", id_1, "vs", id_2, ".smooth.png"))
  smoothScatter(variants$VAF.x, variants$VAF.y, bandwith = 0.1,
                xlab = xl, ylab = yl)
  dev.off()
  
}

plot_all_VAF <- function(samples, mafs, min_depth, min_vaf, genes, effects, out_dir) {
  # Given a vector of samples, plot VAF vs. VAF for all pairs
  
  for (i in 1:length(samples)) {
    
    for (j in 1:length(samples)) {
      
      if (samples[i] == samples[j]) next
      
      maf_1 <- samples[i] %>% read_maf(mafs[i], min_depth)
      maf_2 <- samples[j] %>% read_maf(mafs[j], min_depth)
      
      get_maf_overlap(maf_1, maf_2, min_vaf) %>% 
        plot_VAF_vs_VAF(patient, samples[i], samples[j], genes, effects, out_dir)
    }
  }
}

get_aggregate_maf_df <- function(samples, mafs, min_vaf, min_depth) {
  # Given >=2 MAF paths, return a maf_df which aggregates mutations
  # from all samples
  
  names(mafs) <- samples[seq_along(samples)]
  ag_maf_df <- mafs %>% ldply(read_maf, min_depth, .id = "Sample") %>%                     # "Concatenated" maf (long)
    dplyr::mutate(Sample = paste0(Sample, ".VAF")) %>%                                     # Append VAF
    dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, # Select cols
                  Tumor_Seq_Allele2, Variant_Classification, Sample, VAF) %>% 
    tidyr::spread(Sample, VAF)                                                             # Merge identical variants (wide)
    
  ag_maf_df$maxVAF <- apply(dplyr::select(ag_maf_df, matches(".VAF")), 1, max)
  ag_maf_df %<>% dplyr::filter(maxVAF >= min_vaf) 
  
  return(ag_maf_df[, -which(names(ag_maf_df) == "maxVAF")])

}

plot_vaf_density <- function(patient, samples, mafs, min_vaf, min_depth, out_dir) {
  # Plot VAF density, colour by sample
  
  ag_maf_df <- samples %>% get_aggregate_maf_df(mafs, min_vaf, min_depth) %>% 
    tidyr::gather(Sample, VAF, 8:length(.)) %>%                          # Collapse VAF calls (long)
    tidyr::separate(Sample, into = "sample", sep = "\\.", remove = TRUE) # Remove .VAF suffix from sample col
  
  t_lab <- paste0(patient, "VAF density")
  
  # Plot and save density
  p <- ggplot(ag_maf_df, aes(as.numeric(VAF), colour = sample)) +
    geom_density() + labs(title = t_lab) + xlab("VAF")
  ggsave(p, file = paste0(out_dir, "/", patient, ".density.png"))
  
}

get_private_mutations <- function(patient, samples, mafs, min_vaf, min_depth, genes) {
  # Write private mutations (those with VAF > threshold in
  # one sample but < threshold in all others)
  
  ag_maf_df <- samples %>% get_aggregate_maf_df(mafs, min_vaf, min_depth)
  VAFcols <- names(dplyr::select(ag_maf_df, matches("VAF")))
  
  for (i in 1:length(samples)) {
    
    # Subset to mutations which are above threshold in sample X
    private <- ag_maf_df %>% subset(get(VAFcols[i]) > threshold)    
    
    for (j in 1:(length(VAFcols)-1)) {                
      # Further subset to mutations which are below threshold in all samples but sample X                           
      private <- private %>% subset(get(VAFcols[-i][j]) < threshold)
    }
    
    private %>% write_tsv(file = paste0(out_dir, "/", samples[i], ".private.txt"))
    private %>% dplyr::filter(Hugo_Symbol %in% genes) %>%
      write_tsv(file = paste0(out_dir, "/", samples[i], ".private.in.imp.genes.txt"))
    
  }
}

get_vaf_stats <- function(patient, samples, mafs, min_vaf, min_depth, out_dir) {
  # Get mean/median VAF of all non-private mutations for each sample.
  # A non-private mutation has VAF > min_vaf in all samples
  
  vaf_stats <- samples %>%
    get_aggregate_maf_df(mafs, min_vaf, min_depth) %>%  
    tidyr::gather(Sample, VAF, 8:length(.)) %>%               # Melt VAF column (wide to long format)
    dplyr::filter(VAF > min_vaf) %>%                          # Threshold
    tidyr::separate(Sample, into = "Sample", sep = "\\.",     # Remove .VAF suffix from sample col           
                    extra = "drop", remove = TRUE) %>% 
    dplyr::group_by(Sample) %>% 
    dplyr::summarise(
      Count = n(),
      Mean_VAF = mean(VAF),
      Median_VAF = median(VAF))
  
  vaf_stats %>% write_tsv(file = paste0(out_dir, "/", patient, ".vafstats.tsv"))
    
}