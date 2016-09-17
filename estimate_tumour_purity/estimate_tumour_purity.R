#!/usr/bin/env Rscript
library(MASS)
input_args <- commandArgs(trailingOnly=TRUE)

if (length(input_args) != 1){
  stop("Missing argument: Strelka all.somatic VCF file\n",call.=FALSE)
}

snv_file <- input_args[1]
#snv_file <- "projects/nhl_meta_analysis/all_vcf/MCL_Morin_2015/MCL-13T1_MCL-13N.all.somatic.snvs.vcf"
vcf.colnames <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NORMAL","TUMOR")

somatic.snv.vcf <- read.delim(snv_file,header=F,comment.char='#',col.names=vcf.colnames)

expand_strelka_vcf <- function(vcf){
  
  # Expand INFO column into separate columns for each field
  NT <- character()
  QSS <- numeric()
  QSS_NT <- numeric()
  SGT <- character()
  TQSS <- numeric()
  TQSS_NT <- numeric()
  
  # Expand the normal FORMAT column into separate columns for each field
  N_DP <- numeric()
  N_FDP <- numeric()
  N_SDP <- numeric()
  N_SUBDP <- numeric()
  N_AU_1 <- numeric()
  N_AU_2 <- numeric()
  N_CU_1 <- numeric()
  N_CU_2 <- numeric()
  N_GU_1 <- numeric()
  N_GU_2 <- numeric()
  N_TU_1 <- numeric()
  N_TU_2 <- numeric()
  
  # Expand the tumour FORMAT column into separate columns for each field
  T_DP <- numeric()
  T_FDP <- numeric()
  T_SDP <- numeric()
  T_SUBDP <- numeric()
  T_AU_1 <- numeric()
  T_AU_2 <- numeric()
  T_CU_1 <- numeric()
  T_CU_2 <- numeric()
  T_GU_1 <- numeric()
  T_GU_2 <- numeric()
  T_TU_1 <- numeric()
  T_TU_2 <- numeric()
  
  for(row_n in 1:nrow(vcf)){
    # Order NT,QSS,QSS_NT,SGT,"SOMATIC",TQSS,TQSS_NT
    info.field_values <- unlist(strsplit(as.character(vcf$INFO[row_n]),";"))
    NT <- c(NT, unlist(strsplit(info.field_values[1],"="))[2])
    QSS <- c(QSS, unlist(strsplit(info.field_values[2],"="))[2])
    QSS_NT <- c(QSS_NT, as.numeric(unlist(strsplit(info.field_values[3],"=")))[2])
    SGT <- c(SGT, unlist(strsplit(info.field_values[4],"="))[2])
    TQSS <- c(TQSS, unlist(strsplit(info.field_values[6],"="))[2])
    TQSS_NT <- c(TQSS_NT, unlist(strsplit(info.field_values[7],"="))[2])
    
    # Order N_DP,N_FDP,N_SDP,N_SUBDP,N_AU_1,N_AU_2,N_CU_1,N_CU_2,N_GU_1,N_GU_2,N_TU_1,N_TU_2
    normal.format_values <- unlist(strsplit(as.character(vcf$NORMAL[row_n]),":"))
    N_DP <- c(N_DP, normal.format_values[1])
    N_FDP <- c(N_FDP, normal.format_values[2])
    N_SDP <- c(N_SDP, normal.format_values[3])
    N_SUBDP <- c(N_SUBDP, normal.format_values[4])
    N_AU_1 <- c(N_AU_1, as.numeric(unlist(strsplit(as.character(normal.format_values[5]),","))[1]))
    N_AU_2 <- c(N_AU_2, as.numeric(unlist(strsplit(as.character(normal.format_values[5]),","))[2]))
    N_CU_1 <- c(N_CU_1, as.numeric(unlist(strsplit(as.character(normal.format_values[6]),","))[1]))
    N_CU_2 <- c(N_CU_2, as.numeric(unlist(strsplit(as.character(normal.format_values[6]),","))[2]))
    N_GU_1 <- c(N_GU_1, as.numeric(unlist(strsplit(as.character(normal.format_values[7]),","))[1]))
    N_GU_2 <- c(N_GU_2, as.numeric(unlist(strsplit(as.character(normal.format_values[7]),","))[2]))
    N_TU_1 <- c(N_TU_1, as.numeric(unlist(strsplit(as.character(normal.format_values[8]),","))[1]))
    N_TU_2 <- c(N_TU_2, as.numeric(unlist(strsplit(as.character(normal.format_values[8]),","))[2]))
    
    # Order T_DP,T_FDP,T_SDP,T_SUBDP,T_AU_1,T_AU_2,T_CU_1,T_CU_2,T_GU_1,T_GU_2,T_TU_1,T_TU_2
    tumour.field_values <- unlist(strsplit(as.character(vcf$TUMOR[row_n]),":"))
    T_DP <- c(T_DP, tumour.field_values[1])
    T_FDP <- c(T_FDP, tumour.field_values[2])
    T_SDP <- c(T_SDP, tumour.field_values[3])
    T_SUBDP <- c(T_SUBDP, tumour.field_values[4])
    T_AU_1 <- c(T_AU_1, as.numeric(unlist(strsplit(as.character(tumour.field_values[5]),","))[1]))
    T_AU_2 <- c(T_AU_2, as.numeric(unlist(strsplit(as.character(tumour.field_values[5]),","))[2]))
    T_CU_1 <- c(T_CU_1, as.numeric(unlist(strsplit(as.character(tumour.field_values[6]),","))[1]))
    T_CU_2 <- c(T_CU_2, as.numeric(unlist(strsplit(as.character(tumour.field_values[6]),","))[2]))
    T_GU_1 <- c(T_GU_1, as.numeric(unlist(strsplit(as.character(tumour.field_values[7]),","))[1]))
    T_GU_2 <- c(T_GU_2, as.numeric(unlist(strsplit(as.character(tumour.field_values[7]),","))[2]))
    T_TU_1 <- c(T_TU_1, as.numeric(unlist(strsplit(as.character(tumour.field_values[8]),","))[1]))
    T_TU_2 <- c(T_TU_2, as.numeric(unlist(strsplit(as.character(tumour.field_values[8]),","))[2]))
  }
  
  expanded.vcf <- cbind.data.frame(somatic.snv.vcf,NT,QSS,QSS_NT,SGT,TQSS,TQSS_NT,N_DP,N_FDP,N_SDP,N_SUBDP,N_AU_1,N_AU_2,N_CU_1,N_CU_2,N_GU_1,N_GU_2,N_TU_1,N_TU_2,T_DP,T_FDP,T_SDP,T_SUBDP,T_AU_1,T_AU_2,T_CU_1,T_CU_2,T_GU_1,T_GU_2,T_TU_1,T_TU_2)
  
  return(expanded.vcf)
}

calculate_alt_vafs_in_normal <- function(expanded.vcf){
  # Plot the ALT allele VAF in normal
  # 1. Determine which tier is used for QSS_NT
  #   1.1. Calculate the total read depth in normal depending on the data tier (Skip totals <= 0)
  #     1.1.1. (tier 1) Total = Tier 1 REF count + Tier 2 ALT count
  #     1.1.2. (tier 2) Total = Tier 2 REF count + Tier 2 ALT count
  # 2. Find VAF of ALT allele in normal (Skip ALT columns with two values or if ALT == ".")
  #   2.1. Get correct ALT allele read count depending on the data tier
  #     2.1.1. (tier 1) First number
  #     2.1.2. (tier 2) Second number
  #   2.2. VAF = ALT allele read count / Total
  VAFs <- numeric()
  for(row_n in 1:nrow(expanded.vcf)){
    tqss_nt <- as.numeric(expanded.vcf$TQSS_NT[row_n])
    
    ref_allele <- as.character(expanded.vcf$REF[row_n])
    alt_allele <- as.character(expanded.vcf$ALT[row_n])
    
    if(length(unlist(strsplit(alt_allele,","))) > 1 | alt_allele == "." ){
      next
    }
    
    ref_count <- NULL
    alt_count <- NULL
    if(tqss_nt == 1){
      if(ref_allele == "A"){
        ref_count <- as.numeric(expanded.vcf$N_AU_1[row_n])
      }else if(ref_allele == "C"){
        ref_count <- as.numeric(expanded.vcf$N_CU_1[row_n])
      }else if(ref_allele == "G"){
        ref_count <- as.numeric(expanded.vcf$N_GU_1[row_n])
      }else if(ref_allele == "T"){
        ref_count <- as.numeric(expanded.vcf$N_TU_1[row_n])
      }
      if(alt_allele == "A"){
        alt_count <- as.numeric(expanded.vcf$N_AU_1[row_n])
      }else if(alt_allele == "C"){
        alt_count <- as.numeric(expanded.vcf$N_CU_1[row_n])
      }else if(alt_allele == "G"){
        alt_count <- as.numeric(expanded.vcf$N_GU_1[row_n])
      }else if(alt_allele == "T"){
        alt_count <- as.numeric(expanded.vcf$N_TU_1[row_n])
      }
    }else if(tqss_nt == 2){
      if(ref_allele == "A"){
        ref_count <- as.numeric(expanded.vcf$N_AU_2[row_n])
      }else if(ref_allele == "C"){
        ref_count <- as.numeric(expanded.vcf$N_CU_2[row_n])
      }else if(ref_allele == "G"){
        ref_count <- as.numeric(expanded.vcf$N_GU_2[row_n])
      }else if(ref_allele == "T"){
        ref_count <- as.numeric(expanded.vcf$N_TU_2[row_n])
      }
      if(alt_allele == "A"){
        alt_count <- as.numeric(expanded.vcf$N_AU_2[row_n])
      }else if(alt_allele == "C"){
        alt_count <- as.numeric(expanded.vcf$N_CU_2[row_n])
      }else if(alt_allele == "G"){
        alt_count <- as.numeric(expanded.vcf$N_GU_2[row_n])
      }else if(alt_allele == "T"){
        alt_count <- as.numeric(expanded.vcf$N_TU_2[row_n])
      }
    }
    
    total_depth <- ref_count + alt_count
    
    if(total_depth <= 0){
      VAFs <- c(VAFs, NA)
      next
    }
    
    VAF <- alt_count/total_depth
    VAFs <- c(VAFs, VAF)
    
  }
  return(VAFs)
}

calculate_alt_vafs_in_tumour <- function(expanded.vcf){
  # Plot the ALT allele VAF in normal
  # 1. Determine which tier is used for QSS_NT
  #   1.1. Calculate the total read depth in normal depending on the data tier (Skip totals <= 0)
  #     1.1.1. (tier 1) Total = Tier 1 REF count + Tier 2 ALT count
  #     1.1.2. (tier 2) Total = Tier 2 REF count + Tier 2 ALT count
  # 2. Find VAF of ALT allele in tumour (Skip ALT columns with two values or if ALT == ".")
  #   2.1. Get correct ALT allele read count depending on the data tier
  #     2.1.1. (tier 1) First number
  #     2.1.2. (tier 2) Second number
  #   2.2. VAF = ALT allele read count / Total
  VAFs <- numeric()
  for(row_n in 1:nrow(expanded.vcf)){
    tqss_nt <- as.numeric(expanded.vcf$TQSS_NT[row_n])
    
    ref_allele <- as.character(expanded.vcf$REF[row_n])
    alt_allele <- as.character(expanded.vcf$ALT[row_n])
    
    if(length(unlist(strsplit(alt_allele,","))) > 1 | alt_allele == "." ){
      next
    }
    
    ref_count <- NULL
    alt_count <- NULL
    if(tqss_nt == 1){
      if(ref_allele == "A"){
        ref_count <- as.numeric(expanded.vcf$T_AU_1[row_n])
      }else if(ref_allele == "C"){
        ref_count <- as.numeric(expanded.vcf$T_CU_1[row_n])
      }else if(ref_allele == "G"){
        ref_count <- as.numeric(expanded.vcf$T_GU_1[row_n])
      }else if(ref_allele == "T"){
        ref_count <- as.numeric(expanded.vcf$T_TU_1[row_n])
      }
      if(alt_allele == "A"){
        alt_count <- as.numeric(expanded.vcf$T_AU_1[row_n])
      }else if(alt_allele == "C"){
        alt_count <- as.numeric(expanded.vcf$T_CU_1[row_n])
      }else if(alt_allele == "G"){
        alt_count <- as.numeric(expanded.vcf$T_GU_1[row_n])
      }else if(alt_allele == "T"){
        alt_count <- as.numeric(expanded.vcf$T_TU_1[row_n])
      }
    }else if(tqss_nt == 2){
      if(ref_allele == "A"){
        ref_count <- as.numeric(expanded.vcf$T_AU_2[row_n])
      }else if(ref_allele == "C"){
        ref_count <- as.numeric(expanded.vcf$T_CU_2[row_n])
      }else if(ref_allele == "G"){
        ref_count <- as.numeric(expanded.vcf$T_GU_2[row_n])
      }else if(ref_allele == "T"){
        ref_count <- as.numeric(expanded.vcf$T_TU_2[row_n])
      }
      if(alt_allele == "A"){
        alt_count <- as.numeric(expanded.vcf$T_AU_2[row_n])
      }else if(alt_allele == "C"){
        alt_count <- as.numeric(expanded.vcf$T_CU_2[row_n])
      }else if(alt_allele == "G"){
        alt_count <- as.numeric(expanded.vcf$T_GU_2[row_n])
      }else if(alt_allele == "T"){
        alt_count <- as.numeric(expanded.vcf$T_TU_2[row_n])
      }
    }
    
    total_depth <- ref_count + alt_count
    
    if(total_depth <= 0){
      VAFs <- c(VAFs, NA)
      next
    }
    
    VAF <- alt_count/total_depth
    VAFs <- c(VAFs, VAF)
    
  }
  return(VAFs)
}

somatic.snv.vcf.XInfo.XNormFormat.XTumFormat <- expand_strelka_vcf(somatic.snv.vcf)
somatic.snv.vcf.XInfo.XNormFormat.XTumFormat.alt_filtered1 <- somatic.snv.vcf.XInfo.XNormFormat.XTumFormat[!seq(1,nrow(somatic.snv.vcf.XInfo.XNormFormat.XTumFormat)) %in% grep(",",somatic.snv.vcf.XInfo.XNormFormat.XTumFormat$ALT),]
snv.filtered <- somatic.snv.vcf.XInfo.XNormFormat.XTumFormat.alt_filtered1[somatic.snv.vcf.XInfo.XNormFormat.XTumFormat.alt_filtered1$ALT != ".",]

normal.pass.alt_vaf <- calculate_alt_vafs_in_normal(snv.filtered)
tumour.pass.alt.vaf <- calculate_alt_vafs_in_tumour(snv.filtered)
snv.filtered.alt_vaf <- cbind.data.frame(snv.filtered,vaf_norm=normal.pass.alt_vaf,vaf_tumour=tumour.pass.alt.vaf)

het_put <- snv.filtered.alt_vaf[snv.filtered.alt_vaf$vaf_norm >= 0.2 & snv.filtered.alt_vaf$vaf_norm <= 0.8 & as.numeric(as.character(snv.filtered.alt_vaf$N_DP)) > 20 & as.numeric(as.character(snv.filtered.alt_vaf$T_DP)) > 20,]
#plot(density(het_put$vaf_norm[!is.na(het_put$vaf_norm)]))
# Fit normal distribution on putative SNPs
het_distr <- coef(fitdistr(het_put[!is.na(het_put$vaf_norm),"vaf_norm"],"normal"))
# Get lower and upper threshold to capture 95% of variants; designate these as heterozygous SNPs
het_thresh <- qnorm(p=c(0.025,0.975),mean=het_distr[1],sd=het_distr[2])
het <- snv.filtered.alt_vaf[snv.filtered.alt_vaf$vaf_norm >= het_thresh[1] & snv.filtered.alt_vaf$vaf_norm <= het_thresh[2],]
het.na <- het[!is.na(het$vaf_norm),]
#plot(density(het.na[,"vaf_norm"]))
#plot(density(het.na[,"vaf_tumour"]))
loh.var <- het.na[het.na$vaf_tumour > 0.5,]
#plot(density(loh.var$vaf_tumour))
# Fit normal distribution on LOH events
loh_distr <- coef(fitdistr(loh.var$vaf_tumour,"normal"))
# Use mean as estimate for tumour purity
cat(loh_distr[1])
