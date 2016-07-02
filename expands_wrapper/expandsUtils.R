# Load and process Titan output file
# Header:
# Sample	Chromosome	Start_Position(bp)	End_Position(bp)	Length(bp)
# Median_Ratio	Median_logR	TITAN_state	TITAN_call	Copy_Number	MinorCN
# MajorCN	Clonal_Cluster	Clonal_Frequency
process_titan_seg <- function(seg, include_loh, cn_style) {
  
  print(paste("Loading Titan output file ", seg))
  
  seg1 <- read.csv(seg, stringsAsFactors = FALSE, sep = '\t')
  loh_string = "no_LOH_"
  loh_snv_data <- NULL
  
  if (include_loh > 0) {
  
    if (include_loh == 1) {
      
      neut <- seg1[, "Copy_Number"] == 2
      loh1 <- seg1[, "MajorCN"] == 2
      keep_loh <- neut & loh1
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
    
    sequenza_loh <- seg1[keep_loh, ]
    
    vaf_loh_event <- sequenza_loh[, "Median_Ratio"]  
    loh_snv_data <- matrix(nrow = length(vaf_loh_event), ncol = 7, 
                           dimnames = list(c(), c("chr", "startpos", "endpos", "REF", "ALT", "AF_Tumor", "PN_B")))
    loh_snv_data[, "chr"]      <- as.numeric(sequenza_loh[, "Chromosome"])
    loh_snv_data[, "startpos"] <- as.numeric(sequenza_loh[, 3] + 1)
    loh_snv_data[, "endpos"]   <- as.numeric(sequenza_loh[, 3] + 1)
    loh_snv_data[, "AF_Tumor"] <- as.numeric(sequenza_loh[, "Median_Ratio"])
    loh_snv_data[, "REF"]      <- as.numeric(rep(65, length(vaf_loh_event)))
    loh_snv_data[, "ALT"]      <- as.numeric(rep(45, length(vaf_loh_event)))
    loh_snv_data[, "PN_B"]     <- as.numeric(rep(1, length(vaf_loh_event)))
    
  }
  
  chroms <- as.numeric(seg1[, "Chromosome"])
  
  keep_seg <- seg1[ ,4] - seg1[, 3] > 1
  
  keep_chrom <- !is.na(chroms < 23)
  
  keep_both <- keep_chrom & keep_seg
  seg2 <- matrix(nrow = dim(seg1[keep_both, ][1]), ncol = 4)
  
  colnames(seg2) <- c("chr", "startpos", "endpos", "CN_Estimate")
  
  seg2[, "chr"] <- as.numeric(seg1[keep_both, "Chromosome"])
  
  seg2[, "startpos"] <- as.numeric(seg1[keep_both, 3])
  seg2[, "endpos"]   <- as.numeric(seg1[keep_both, 4])
  
  if (cn_style == 1) {
    seg2[, 4] <- as.numeric(seg1[keep_both, "Copy_Number"])
  }
  else if (cn_style == 2) {
    seg2[ ,4] <- 2*2^seg1[keep_both, "Median_logR"]
  }	
  
  output <- list(seg2, loh_snv_data)
  return(output)
  
}

# Expecteds Sequenza-style seg file
# Header:
# chromosome	start.pos	 end.pos	Bf	N.BAF	 sd.BAF	depth.ratio	 
# N.ratio	 sd.ratio	 CNt 	A	  B	 LPP
process_sequenza_seg <- function(seg, include_loh, cn_style) {
  
  print(paste("Loading Sequenza seg file ", seg))
  seg1 <- read.csv(seg, stringsAsFactors = FALSE, sep = '\t')
  # note that the CNt value is being used directly by Expands.
  # For Titan, one would have to convert the log ratio to absolute copy number
  
  loh_string = "no_LOH_"
  loh_snv_data <- NULL
  
  if (include_loh > 0) {
    
    if (include_loh == 1) {
      
      neut <- seg1[, "CNt"] == 2
      loh1 <- seg1[, "A"] == 2
      keep_loh <- neut & loh1
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
  
  chroms_a <- seg1[, 2]
  chroms <- as.numeric(chroms_a)
  
  keep_seg <- seg1[, "end.pos"] - seg1[, "start.pos"] > 1
  
  keep_chrom <- !is.na(chroms < 23)
  
  keep_both <- keep_chrom & keep_seg
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

# Expects OncoSNP .cnvs file augmented with Log R Ratio and BAF using oncosnputils
# Columns: chr	start	end	copyNum	loh	rank	logLik	numMarkers	normFrac	state	ploidyNum
# majorCopyNumber minorCopyNumber	LRR	LRRShifted	BAF	numProbes	numSnpProbes
process_oncosnp_seg <- function(seg, include_loh, cn_style) {
  
  print(paste("loading OncoSNP ", seg))
  seg1 <- read.csv(seg, stringsAsFactors = FALSE, sep = "\t")
  chroms_a  <- seg1[, 1]
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
  
  snv_data <- matrix(nrow = dim(maf_keep[1]), ncol = 7,
                     dimnames = list(c(), c("chr", "startpos", "endpos", "REF", "ALT", "AF_Tumor", "PN_B")))
  
  # Remove chr prefix
  snv_data[, "chr"] <- as.numeric(sub("^chr", "", maf_data[, "Chromosome"]))
  
  # Load start and end position into matrix
  snv_data[, "startpos"] <- as.numeric(maf_keep[, "Start_Position"])
  snv_data[, "endpos"]   <- as.numeric(maf_keep[, "End_Position"])
  snv_data[, "REF"]      <- sapply(maf_keep[, "Reference_Allele"], utf8ToInt)
  snv_data[, "ALT"]      <- sapply(maf_keep[, "Tumor_Seq_Allele2"], utf8ToInt)
  snv_data[, "AF_Tumor"] <- maf_keep[, "t_alt_count"] / (maf_keep[, "t_alt_count"] + maf_keep[, "t_ref_count"])
  
  # set flag to somatic for all SNVs
  snv_data[, "PN_B"] <- 0
  
  outputs <- list(snv_data, maf_keep)
  return(outputs)
  
}

# Assign copy numbers in cbs to mutations in dm
assign_states_to_mutation <- function(dm, cbs, cols) {
  print("Assigning copy number to mutations for PyClone... ")
  
  for (k in 1:nrow(cbs)){
    
    # Get all mutations in current segment
    idx <- which(as.numeric(dm[, "chr"]) == cbs[k, "chr"] & as.numeric(dm[, "startpos"]) >= cbs[k, "startpos"] & as.numeric(dm[, "startpos"]) <= cbs[k, "endpos"])
    
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
    idx <- which(as.numeric(dm[k, "chr"]) == cbs[, "chr"] & as.numeric(dm[k, "startpos"]) >= cbs[, "startpos"] & as.numeric(dm[k, "endpos"]) <= cbs[, "endpos"])
    
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
  
  dm <- dm[, colnames(dm) != "segmentLength"]
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
  
  # Remove chr prefix
  py_snv_data[, "chr"] <- as.numeric(sub("^chr", "", maf_keep[,"Chromosome"]))
  
  # load start and end position into matrix
  py_snv_data[, "startpos"]   <- as.numeric(maf_keep[, "Start_Position"])
  py_snv_data[, "endpos"]     <- as.numeric(maf_keep[, "End_Position"])
  
  py_snv_data[, "ref_counts"] <- maf_keep[, "t_ref_count"]
  py_snv_data[, "var_counts"] <- maf_keep[, "t_alt_count"]
  
  # seg rename columns
  if(input_mode == "S") {
    colnames(seg1) <- c("chr", "startpos", "endpos", "Bf", "N.BAF", "sd.BAF", "depth.ratio",
                        "N.ratio", "sd.ratio", "CNt", "major_cn", "minor_cn", "segmentLength")
    
  } else if (input_mode == "O") {
    colnames(seg1) <- c("chr", "startpos", "endpos", "CN", "loh", "rank", "logLik", "numMarkers", 
                        "normFrac", "state", "ploidyNum", "major_cn", "minor_cn", "log.ratio",
                        "log.ratio.shifted", "BAF", "numProbes", "numSnpProbes")
  }
  
  seg1[, "normal_cn"] <- 2
  seg1[, "segmentLength"] <- seg1[, "endpos"] - seg1[, "startpos"]
  
  # if oncosnp, assign states by rank, otherwise, use general assign_states function
  if (input_mode == "O") {
    py_snv_data_assigned <- assign_states_to_mutation_by_rank(py_snv_data, seg1, c("minor_cn", "major_cn"))
  } else {
    py_snv_data_assigned <- assign_states_to_mutation(py_snv_data, seg1, c("minor_cn", "major_cn"))
  }
  
  py_snv_data_assigned[, "gene"]        <- maf_keep[, "Hugo_Symbol"]
  py_snv_data_assigned[, "mutation_id"] <- paste(py_snv_data_assigned[, "gene"],
                                                 py_snv_data_assigned[, "startpos"], sep = "_")
  
  keepers <- !is.na(py_snv_data_assigned[, "normal_cn"])
  
  # preserve structure even if only 1 variant returned
  return(py_snv_data_assigned[keepers, , drop = FALSE])
  
}
