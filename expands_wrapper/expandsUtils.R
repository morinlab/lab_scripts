# Load and process Titan output file
# Sample	Chromosome	Start_Position(bp)	End_Position(bp)	Length(bp)
# Median_Ratio	Median_logR	TITAN_state	TITAN_call	Copy_Number	MinorCN
# MajorCN	Clonal_Cluster	Clonal_Frequency
process_titan_seg <- function(seg, include_loh, cn_style) {
  
  print(paste("Loading Titan output file ", seg))
  
  seg1 <- read.csv(seg, stringsAsFactors = FALSE, sep = '\t')
  
  # Remove chr prefix if it exists
  test_chr <- seg1[1, "Chromosome"]
  if (grepl("chr", test_chr)) {
    # 'chr' prefix exists
    seg1[, "Chromosome"] <- sub("^chr", "", seg1[, "Chromosome"])
  }
  chroms <- as.numeric(seg1[, "Chromosome"])
  keep_chrom <- !is.na(chroms < 23)
  
  loh_string = "no_LOH_"
  loh_snv_data <- NULL
  
  if (include_loh > 0) {
  
    if (include_loh == 1) {
      
      neut <- seg1[, "Copy_Number"] == 2
      loh1 <- seg1[, "MajorCN"] == 2
      keep_loh <- neut & loh1 & keep_chrom
      loh_string <- "neut_LOH_"
      
    } else if (include_loh == 2) {
    
      # warning: These other two options should not be used. I don't think they work well. Please stick with 0 or 1
      del <- seg1[, "Copy_Number"] == 1
      loh1 <-  seg1[, "MajorCN"] + seg1[, "MinorCN"] == 1
      keep_loh <- del & loh1
      loh_string <- "del_LOH_"
      
    } else if (include_loh == 3) {
      
      loh_string <- "any_LOH_"
      not_amp <- seg1[,"Copy_Number"] <= 2
      loh1 <- seg1[,"MinorCN"]==0
      keep_loh <- not_amp & loh1
      
    }
    
    is_loh <- seg1[keep_loh, ]
    
    vaf_loh_event <- is_loh[, "Median_Ratio"]  
    loh_snv_data <- matrix(nrow = length(vaf_loh_event), ncol = 7, 
                           dimnames = list(c(), c("chr", "startpos", "endpos", "REF", "ALT", "AF_Tumor", "PN_B")))
    loh_snv_data[, "chr"]      <- as.numeric(is_loh[, "Chromosome"])
    loh_snv_data[, "startpos"] <- as.numeric(is_loh[, 3] + 1)
    loh_snv_data[, "endpos"]   <- as.numeric(is_loh[, 3] + 1)
    loh_snv_data[, "AF_Tumor"] <- as.numeric(is_loh[, "Median_Ratio"])
    loh_snv_data[, "REF"]      <- as.numeric(rep(65, length(vaf_loh_event)))
    loh_snv_data[, "ALT"]      <- as.numeric(rep(45, length(vaf_loh_event)))
    loh_snv_data[, "PN_B"]     <- as.numeric(rep(1, length(vaf_loh_event)))
  }
  
  # Segment length > 1 and autosomes only
  keep_seg <- seg1[ ,4] - seg1[, 3] > 1
  keep_both <- keep_chrom & keep_seg
  
  # Make input CNV matrix for EXPANDS
  seg2 <- matrix(nrow = dim(seg1[keep_both, ][1]), ncol = 4)
  colnames(seg2) <- c("chr", "startpos", "endpos", "CN_Estimate")
  
  seg2[, "chr"] <- as.numeric(seg1[keep_both, "Chromosome"])
  seg2[, "startpos"] <- as.numeric(seg1[keep_both, 3])
  seg2[, "endpos"]   <- as.numeric(seg1[keep_both, 4])
  
  if (cn_style == 1) {
    seg2[, 4] <- as.numeric(seg1[keep_both, "Copy_Number"])
  } else if (cn_style == 2) {
    seg2[ ,4] <- 2*2^seg1[keep_both, "Median_logR"]
  }	
  
  output <- list(seg2, loh_snv_data)
  return(output)
  
}

# Load and process Sequenza-style seg file
# chromosome	start.pos	 end.pos	Bf	N.BAF	 sd.BAF	depth.ratio	 
# N.ratio	 sd.ratio	 CNt 	A	  B	 LPP
process_sequenza_seg <- function(seg, include_loh, cn_style) {
  
  print(paste("Loading Sequenza seg file ", seg))
  seg1 <- read.csv(seg, stringsAsFactors = FALSE, sep = '\t')
  
  # Remove chr prefix if it exists
  test_chr <- seg1[1, "chromosome"]
  if (grepl("chr", test_chr)) {
    # 'chr' prefix exists
    seg1[, "chromosome"] <- sub("^chr", "", seg1[, "chromosome"])
  }
  chroms_a <- seg1[, "chromosome"]
  chroms <- as.numeric(chroms_a)
  keep_chrom <- !is.na(chroms < 23)
  
  loh_string = "no_LOH_"
  loh_snv_data <- NULL
  
  if (include_loh > 0) {
    
    if (include_loh == 1) {
      
      neut <- seg1[, "CNt"] == 2
      loh1 <- seg1[, "A"] == 2
      keep_loh <- neut & loh1 & keep_chrom
      loh_string <- "neut_LOH_"
      
    } else if (include_loh == 2) {
      
      #don't use these other two options!
      del <- seg1[, "CNt"] == 1
      loh1 <- seg1[, "A"] + seg1[, "B"] == 1
      keep_loh <- del & loh1
      loh_string <- "del_LOH_"
      
    } else if (include_loh == 3) {
      
      loh_string <- "any_LOH"
      not_amp <- seg1[, "CNt"] <= 2
      loh1 <- seg1[, "B"] == 0
      keep_loh <- not_amp & loh1
    }
    
    sequenza_loh <- seg1[keep_loh, ]
    
    vaf_loh_event <- 1- sequenza_loh[, "Bf"]
    # the Expands method expects VAFs that were increased relative to the normal
    # but Sequenza reports BAF relative to lower Reference_Allele
    
    loh_snv_data <- matrix(nrow = length(vaf_loh_event), ncol = 7,
                        dimnames = list(c(), c("chr", "startpos", "endpos", "REF", "ALT", "AF_Tumor", "PN_B")))
    loh_snv_data[, "chr"]      <- as.numeric(sequenza_loh[, "chromosome"])
    loh_snv_data[, "startpos"] <- as.numeric(sequenza_loh[, "start.pos"] + 1)
    loh_snv_data[, "endpos"]   <- as.numeric(sequenza_loh[, "start.pos"] + 1)
    loh_snv_data[, "AF_Tumor"] <- as.numeric(vaf_loh_event)
    loh_snv_data[, "REF"]      <- as.numeric(rep(65, length(vaf_loh_event)))
    loh_snv_data[, "ALT"]      <- as.numeric(rep(45, length(vaf_loh_event)))
    loh_snv_data[, "PN_B"]     <- as.numeric(rep(1, length(vaf_loh_event)))
    # loh_snv_data[,"PN_B"] = as.numeric(rep(0,length(vaf_loh_event)))
    # don't set as LOH because it works better if you don't. Unsure why

  }
  
  keep_seg <- seg1[, "end.pos"] - seg1[, "start.pos"] > 1
  keep_both <- keep_chrom & keep_seg
  
  # Make input CNV matrix for EXPANDS
  seg2 <- matrix(nrow = dim(seg1[keep_both, ][1]), ncol = 4)
  colnames(seg2) <- c("chr",  "startpos","endpos","CN_Estimate")
  
  seg2[, "chr"]      <- as.numeric(seg1[keep_both, "chromosome"])
  seg2[, "startpos"] <- as.numeric(seg1[keep_both, "start.pos"])
  seg2[, "endpos"]   <- as.numeric(seg1[keep_both, "end.pos"])
  
  if(cn_style == 1) {
    seg2[, 4] <- as.numeric(seg1[keep_both, "CNt"])
  } else if (cn_style == 2) {
    logratio <- log2(as.numeric(seg1[keep_both, "depth.ratio"]))
    seg2[, 4] <- 2*2^logratio
  }
  
  output <- list(seg2, loh_snv_data)
  return(output)

}

# untested
process_igv_seg <- function(seg, include_loh, cn_style) {
  
  print(paste("Loading IGV-friendly seg file ", seg))
  # LOH information needs to be included in these files and BAF
  # example format: sample	chr	start	end	LOH_flag	BAF	median
  
  seg1 <- read.csv(seg, stringsAsFactors = FALSE, sep = '\t')
  # convert last column from log ratio to absolute CN
  logratios <- seg1[, length(seg1[1, ])]
  abs <- 2*2^logratios
  chroms_a <- seg1[,2]
  chroms <- as.numeric(chroms_a)
  keep_chrom <- !is.na(chroms<23)
  
  #loh_string = "no_LOH_"
  loh_snv_data <- NULL
  loh_status = seg1[, "LOH_flag"]
  
  if(include_loh > 0){
    #Warning: This is the most experimental mode of running currently. 
    
    sequenza_loh <- seg1[loh_status == 1, ]
    vaf_loh_event <- 1- sequenza_loh[, "BAF"]
    # the Expands method expects VAFs that were increased relative to the normal but
    # Sequenza reports BAF relative to lower Reference_Allel
    
    loh_snv_data <- matrix(nrow = length(vaf_loh_event), ncol = 7,
                           dimnames = list(c(),c("chr", "startpos", "endpos", "REF", "ALT", "AF_Tumor", "PN_B")))
    loh_snv_data[, "chr"]      <- as.numeric(sequenza_loh[, "chromosome"])
    loh_snv_data[, "startpos"] <- as.numeric(sequenza_loh[, "start"] + 1)
    loh_snv_data[, "endpos"]   <- as.numeric(sequenza_loh[, "end"] + 1)
    loh_snv_data[, "AF_Tumor"] <- as.numeric(vaf_loh_event)
    loh_snv_data[, "REF"]      <- as.numeric(rep(65, length(vaf_loh_event)))
    loh_snv_data[, "ALT"]      <- as.numeric(rep(45, length(vaf_loh_event)))
    loh_snv_data[, "PN_B"]     <- as.numeric(rep(1, length(vaf_loh_event)))
  
  }
  
  seg2 <- matrix(nrow = length(seg1[keep_chrom, 1]), ncol = 4)
  
  #seg2=matrix(nrow=length(seg1[keep_chrom,1]),ncol=4)
  colnames(seg2) <- c("chr", "startpos", "endpos", "CN_Estimate")
  
  seg2[, "chr"]      <- as.numeric(seg1[keep_chrom, 2])
  seg2[, "startpos"] <- as.numeric(seg1[keep_chrom, 3])
  seg2[, "endpos"]   <- as.numeric(seg1[keep_chrom, 4])
  seg2[, 4]          <- abs[keep_chrom]
  
  output <- list(seg2, loh_snv_data)
  return(output)
}

# Load and process OncoSNP .cnvs file augmented with Log R Ratio and BAF using oncosnputils
# chr	start	end	copyNum	loh	rank	logLik	numMarkers	normFrac	state	ploidyNum
# majorCopyNumber minorCopyNumber	LRR	LRRShifted	BAF	numProbes	numSnpProbes
process_oncosnp_seg <- function(seg, include_loh, cn_style) {
  
  print(paste("loading OncoSNP ", seg))
  seg1 <- read.csv(seg, stringsAsFactors = FALSE, sep = "\t")
  chroms_a  <- seg1[, "chr"]
  chroms <- as.numeric(chroms_a)
  
  keep_seg <- seg1[, "end"] - seg1[, "start"] > 1
  keep_chrom  <- !is.na(chroms < 23)
  
  keep_both <- keep_chrom & keep_seg
  seg2 <- matrix(nrow = dim(seg1[keep_both, ][1]), ncol = 4)
  
  colnames(seg2) <- c("chr", "startpos", "endpos", "CN_Estimate")
  
  seg2[, "chr"]      <- as.numeric(seg1[keep_both, "chr"])
  seg2[, "startpos"] <- as.numeric(seg1[keep_both, "start"])
  seg2[, "endpos"]   <- as.numeric(seg1[keep_both, "end"])
  
  if (cn_style == 1) {
    
    seg2[, 4] <- as.numeric(seg1[keep_both, "copyNum"])
    
  } else if (cn_style == 2) {
    
    logratio <- as.numeric(seg1[keep_both, "LRR"])
    seg2[, 4] <- 2*2^logratio
    
  }

  loh_snv_data <- NULL
  
  output <- list(seg2, loh_snv_data)
  return(output)
  
}

process_maf <- function(maf) {
  
  print(paste("Loading MAF ", maf))
  maf_data <- read.csv(maf, sep = '\t', stringsAsFactors = FALSE )
  
  maf_keep_cols <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2",
                    "t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count")
  
  # this method works with Indels but is a bit hacky because a fake ref/alt value is used.
  # It doesn't really matter as long as you ignore those values in the output. Expands ignores them AFAIK
  maf_keep <- maf_data[, maf_keep_cols]
  maf_keep[, "Reference_Allele"]  <- sapply(maf_keep[, "Reference_Allele"], function(x) substr(x,1,1))
  maf_keep[, "Tumor_Seq_Allele2"] <- sapply(maf_keep[, "Tumor_Seq_Allele2"], function(x) substr(x,1,1))
  
  # Remove chr prefix if it exists
  test_chr <- maf_keep[1, "Chromosome"]
  if (grepl("chr", test_chr)) {
    # 'chr' prefix exists
    maf_keep[, "Chromosome"] <- sub("^chr", "", maf_keep[, "Chromosome"])
  }
  
  # Keep mutations on autosomes only
  chroms_a  <- maf_keep[, "Chromosome"]
  chroms <- as.numeric(chroms_a)
  keep_chrom  <- !is.na(chroms < 23)
  
  snv_data <- matrix(nrow = dim(maf_keep[keep_chrom, ][1]), ncol = 7,
                     dimnames = list(c(), c("chr", "startpos", "endpos", "REF", "ALT", "AF_Tumor", "PN_B")))
  
  # Load start and end position into matrix
  snv_data[, "chr"]      <- as.numeric(maf_keep[keep_chrom, "Chromosome"])
  snv_data[, "startpos"] <- as.numeric(maf_keep[keep_chrom, "Start_Position"])
  snv_data[, "endpos"]   <- as.numeric(maf_keep[keep_chrom, "End_Position"])
  snv_data[, "REF"]      <- sapply(maf_keep[keep_chrom, "Reference_Allele"], utf8ToInt)
  snv_data[, "ALT"]      <- sapply(maf_keep[keep_chrom, "Tumor_Seq_Allele2"], utf8ToInt)
  snv_data[, "AF_Tumor"] <- maf_keep[keep_chrom, "t_alt_count"] / (maf_keep[keep_chrom, "t_alt_count"] + maf_keep[keep_chrom, "t_ref_count"])
  
  # Set flag to somatic for all SNVs
  snv_data[, "PN_B"] <- 0
  
  # Exclude mitochondrial chromosomes from maf_keep
  not_chr_m <- !(chroms_a == "M")
  
  outputs <- list(snv_data, maf_keep[not_chr_m, ])
  return(outputs)
  
}

# Assign copy numbers in cbs to mutations in dm
assign_states_to_mutation <- function(dm, cbs, cols) {
  print("Assigning copy number to mutations for PyClone... ")
  
  for (k in 1:nrow(cbs)){
    
    # Get all mutations in current segment
    idx <- which(dm[, "chr"] == cbs[k, "chr"] & as.numeric(dm[, "startpos"]) >= cbs[k, "startpos"] & as.numeric(dm[, "startpos"]) <= cbs[k, "endpos"])
    
    if (length(idx) == 0){
      next
    }
    
    # the only bit of code that relies on the matlab dependency
    dm[idx, cols] <- repmat(as.numeric(cbs[k, cols]), length(idx), 1)
    dm[idx, "normal_cn"] <- rep(2, length(idx))
    
  }
  
  dm <- dm[, colnames(dm) != "segmentLength"]
  print("... Done.")
  
  return(dm)
}

# Assign copy numbers in cbs to mutations in dm
# for OncoSNP output (use max rank)
assign_states_to_mutation_by_rank <- function(dm, cbs, cols) {
  print("Assigning copy number to mutations for PyClone using rank... ")
  
  for (k in 1:nrow(dm)){
    
    best_rank <- 0
    idx <- which(dm[k, "chr"] == cbs[, "chr"] & as.numeric(dm[k, "startpos"]) >= cbs[, "startpos"] & as.numeric(dm[k, "endpos"]) <= cbs[, "endpos"])
    
    if (length(idx) == 0){
      next
    }
    
    # if the mutations is found in multiple segments, assign copy number
    # from the segment with highest rank
    for (j in 1:length(idx)) {
      
        rank <- cbs[j, "rank"]
        
        if (rank > best_rank) {
          top_ranked_idx <- idx[j]
          best_rank <- rank
        }

    }
  
    dm[k, cols] <- repmat(as.numeric(cbs[top_ranked_idx, cols]), 1)
    dm[k, "normal_cn"] <- rep(2, 1)
    
  }
  
  dm <- dm[, colnames(dm) != "segmentLength", drop = FALSE]
  print("... Done.")
  
  return(dm)
}


generate_pyclone_input <- function(seg, maf_keep, input_mode) {

  seg1 <- read.csv(seg, stringsAsFactors = FALSE, sep = "\t")

  # Convert to matrix
  maf_keep <- do.call(rbind, maf_keep)

  py_snv_data <- matrix(nrow = dim(maf_keep[1]), ncol = 10,
                        dimnames = list(c(), c("gene", "chr", "startpos", "endpos", "mutation_id",
                                               "ref_counts", "var_counts", "normal_cn", "major_cn", "minor_cn")))
  
  # Load start and end position into matrix
  py_snv_data[, "startpos"]   <- as.numeric(maf_keep[, "Start_Position"])
  py_snv_data[, "endpos"]     <- as.numeric(maf_keep[, "End_Position"])
  
  # Remove chr prefix
  py_snv_data[, "chr"] <- sub("^chr", "", maf_keep[,"Chromosome"])
  
  py_snv_data[, "ref_counts"] <- maf_keep[, "t_ref_count"]
  py_snv_data[, "var_counts"] <- maf_keep[, "t_alt_count"]
  
  # Rename seg columns
  if(input_mode == "S") {
    colnames(seg1) <- c("chr", "startpos", "endpos", "Bf", "N.BAF", "sd.BAF", "depth.ratio",
                        "N.ratio", "sd.ratio", "CNt", "major_cn", "minor_cn", "segmentLength")
    
  } else if (input_mode == "O") {
    colnames(seg1) <- c("chr", "startpos", "endpos", "CN", "loh", "rank", "logLik", "numMarkers", 
                        "normFrac", "state", "ploidyNum", "major_cn", "minor_cn", "log.ratio",
                        "log.ratio.shifted", "BAF", "numProbes", "numSnpProbes")
  } else if (input_mode == "T") {
    colnames(seg1) <- c("sample", "chr", "startpos", "endpos", "segmentLength", "median.ratio", 
                        "median.log.ratio", "state", "call", "CN", "minor_cn", "major_cn", 
                        "clonal_cluster", "clonal_frequency")
  }
  
  seg1[, "normal_cn"] <- 2
  seg1[, "segmentLength"] <- seg1[, "endpos"] - seg1[, "startpos"]
  
  # If oncosnp, assign states by rank, otherwise, use general assign_states function
  if (input_mode == "O") {
    py_snv_data_assigned <- assign_states_to_mutation_by_rank(py_snv_data, seg1, c("minor_cn", "major_cn"))
  } else {
    py_snv_data_assigned <- assign_states_to_mutation(py_snv_data, seg1, c("minor_cn", "major_cn"))
  }
  
  py_snv_data_assigned[, "gene"]        <- maf_keep[, "Hugo_Symbol"]
  py_snv_data_assigned[, "mutation_id"] <- paste(py_snv_data_assigned[, "gene"],
                                                 py_snv_data_assigned[, "startpos"], sep = "_")
  
  # If no seg data, assume normal CN
  normal_cn_profile <- matrix(nrow = 1, ncol = 3, dimnames = list(c(), c("normal_cn", "minor_cn", "major_cn")))
  normal_cn_profile[, "normal_cn"] <- 2
  normal_cn_profile[, "minor_cn"] <- 0
  normal_cn_profile[, "major_cn"] <- 2
  
  no_cn_data <- is.na(py_snv_data_assigned[, "normal_cn"])
  
  # If any mutations have no segment data, assign normal CN
  if (length(py_snv_data_assigned[no_cn_data, 1]) != 0) {
    
    cols2 = c("normal_cn", "minor_cn", "major_cn")
    py_snv_data_assigned[no_cn_data, cols2] <- repmat(normal_cn_profile[1, cols2], length(py_snv_data_assigned[no_cn_data, 1]), 1)
    
  }
  
  complete <- !is.na(py_snv_data_assigned[, "major_cn"]) & !is.na(py_snv_data_assigned[, "minor_cn"])
  return(py_snv_data_assigned[complete, , drop = FALSE])
  
}

load_plot_libs <- function() {
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(ggrepel))
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(grid))
  suppressPackageStartupMessages(require(gtable))
  suppressPackageStartupMessages(require(dplyr))
}

tidy_dm <- function(dm, orderBy) {
  # Wrangle dm into a tidy data frame for plotting
  
  orderBy <- ifelse(orderBy == "conf", "desc(SP_conf)", orderBy)
  
  var_df <- data.frame(dm) %>%                      # Convert dm to data frame
    filter(!is.na(SP)) %>%                          # Remove variants with no SP assigned
    mutate(SP_conf = X.maxP - min(X.maxP,           # Normalize SP assignment confidence
                                  na.rm = TRUE) + 1) %>% 
    mutate(SP_conf = 50 * SP_conf / max(SP_conf, na.rm = TRUE)) %>%
    arrange_("desc(SP)", orderBy, "startpos") %>%   # Re-order by SP then chr then pos
    mutate(idx = rownames(.))
  
  print("test1")
  return(var_df)
}

adjust_tumor_af <- function(var_df, rawAF) {
  # Given a var_df returned from tidy_dm(), :
  # if plotting rawAF, do no adjustments
  # if plotting adjusted AF, execute adjustment
  
  if (!rawAF) { # Adjust tumor AF
    iEq2 <- which(var_df$PM_B == var_df$PN_B)
    iEq3 <- setdiff(1:nrow(var_df), iEq2)
    
    if (!isempty(iEq3)) {   # Use Equation 3
      var_df[iEq3, "AF_Tumor_Adjusted"] <- (var_df[iEq3, "AF_Tumor"] * var_df[iEq3, "CN_Estimate"] - var_df[iEq3, "PN_B"]) /
        (var_df[iEq3, "PM_B"] - var_df[iEq3, "PN_B"])
      var_df[iEq3, "AF_Tumor_Adjusted"] <- var_df[iEq3, "AF_Tumor_Adjusted"] * (var_df[iEq3, "PM_B"] / var_df[iEq3, "PM"])
    }
    
    if (!isempty(iEq2)) {   # Use Equation 2
      var_df[iEq2, "AF_Tumor_Adjusted"] <- (var_df[iEq2, "CN_Estimate"] - 2) / (var_df[iEq2, "PM"] - 2)
    }
    
  } else {     # Use raw tumor AF
    var_df %<>% mutate(AF_Tumor_Adjusted = AF_Tumor)
  }
  
  print("test2")
  return(var_df)
}

preserve_var_df <- function(var_df) {
  # Turn observations into factors etc
  # to preserve ordering when plotting
  
  var_df$idx <- as.numeric(var_df$idx)
  var_df$idx <- factor(var_df$idx, levels = order(var_df$idx))
  var_df$CN_Estimate <- factor(var_df$CN_Estimate)
  var_df$PN_B <- factor(as.numeric(var_df$PN_B))
  
  print("test3")
  return(var_df)
}

get_data_to_label <- function(var_df, maf, genes, effects) {
  
  maf_df <- read.delim(maf, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  test_chr <- maf_df[1, "Chromosome"]
  if (grepl("chr", test_chr)) {
    # 'chr' prefix exists
    maf_df$Chromosome <- sub("^chr", "", maf_df$Chromosome) %>% as.numeric()
  }
  
  gene_df <- maf_df[, c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification")]
  
  label_df <- left_join(var_df, gene_df,
                         by = c("chr" = "Chromosome", "startpos" = "Start_Position", "endpos" = "End_Position")) %>% 
    dplyr::rename(Gene = Hugo_Symbol) %>% 
    dplyr::rename(Effect = Variant_Classification) %>% 
    dplyr::filter(Gene %in% genes) %>% 
    dplyr::filter(Effect %in% effects)
    
  return(label_df)
}

plot_expands_SPs <- function(dm, sampleID, maf, orderBy = "conf", rawAF = FALSE, genes = NULL, effects) {
  load_plot_libs()
  
  # Get data to plot
  var_df <- tidy_dm(dm, orderBy) %>% 
    adjust_tumor_af(rawAF) %>% 
    preserve_var_df()
  
  # Get data to label
  if (!is.null(genes)) {
    label_df <- get_data_to_label(var_df, maf, genes, effects)
  }
  


  # Use EXPANDS palette for colouring by CN
  cn_palette <- c("1" = "#A6CEE3", "2" = "#1F78B4", "3" = "#B2DF8A", "4" = "#33A02C",
                  "5" = "#FB9A99", "6" = "#E31A1C", "7" = "#FDBF6F", "8" = "#FF7F00")
  yl <- ifelse(rawAF, "Adjusted Allele Frequency", "Allele Frequency")

  # Allele frequency/SP plot
  af_plot <- ggplot(var_df) +
    geom_point(aes(x = idx, y = SP, fill = SP_conf), shape = 22, colour = NA) +
    scale_fill_gradient(low = "#E6E6E6", high = "#4D4D4D", guide = FALSE) +
    geom_point(alpha = 0.8,
               aes(x = idx, y = AF_Tumor_Adjusted, colour = CN_Estimate, shape = PN_B)) +
    scale_shape_discrete(name = "Type", labels = c("SNV", "LOH")) +
    scale_colour_manual(values = cn_palette, name = "Copy\nNumber") +
    ylab(yl) + xlab(NULL) + 
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.line.y = element_line(size = 0.3),
          axis.line.x = element_line(size = 0.3)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) #+ ggtitle(sampleID)

  print("finished af_plot")

  # Add annotation if available
  if (!is.null(genes)) {
    af_plot <- af_plot + geom_text_repel(data = label_df,
                    aes(x = idx, y = AF_Tumor_Adjusted, label = Gene),
                     nudge_x = 200, nudge_y = 0.3, size = 3.5)
  }

  # Copy number/SP plot
  cn_plot <- ggplot(var_df, aes(x = idx, y = PM_cnv)) +
    geom_point(aes(colour = CN_Estimate), alpha = 0.8) +
    scale_colour_manual(values = cn_palette, guide = FALSE) +
    ylab("Total ploidy at locus") + xlab("Mutation") +
    theme(axis.line.x = element_line(size = 0.3),
          axis.line.y = element_line(size = 0.3),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_y_continuous(breaks = seq(1, 5, by = 1), limits = c(1, 5)) +
    scale_x_discrete(breaks = seq(0, nrow(var_df), by = 200)) +
    coord_fixed(ratio = 1000/10)

  print("finished cn_plot")

  # Aggregate plots
  af_g <- ggplotGrob(af_plot)
  cn_g <- ggplotGrob(cn_plot)

  # Add a column of legend width to cn grob
  legend_width <- af_g$widths[5]
  cn_g <- gtable_add_cols(cn_g, legend_width, 4)

  # Add title
  #title <- textGrob(sampleID, gp = gpar(fontsize = 10))
  
  # Bind and arrange the two grobs
  g <- rbind(af_g, cn_g, size = "first")
  g$widths <- unit.pmax(af_g$widths, cn_g$widths)
  #g <- gtable_add_rows(title)
  
  grid.newpage()
  grid.draw(g)

}
