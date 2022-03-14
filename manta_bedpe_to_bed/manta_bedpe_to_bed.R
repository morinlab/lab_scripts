#!/usr/bin/env Rscript
# ----
# Input: Manta BEDPE file
# Example command:
#   Rscript manta_bedpe_to_bed.R hg19 input.bedpe output.pass.bedpe output.pass.bed output.ig_myc.bed
# ----


# Load Packages -------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})


# Argument Parsing ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

genome_build <- args[1]

paths <- list()
paths$input_bedpe     <- args[2]
paths$output_bedpe    <- args[3]
paths$output_bed      <- args[4]
paths$output_bed_myc  <- args[5]


# Set reference based on genome build -------------------------------------

if (genome_build == "hg38") {
  myc_exon_3_midpoint <- 127740675
  loci <- list(
    IGK = list(chrom = "chr2", start = 87999518, end = 90599757),
    MYC = list(chrom = "chr8", start = 126393182, end = 130762146),
    IGH = list(chrom = "chr14", start = 104589639, end = 107810399),
    IGL = list(chrom = "chr22", start = 21031465, end = 23905532))
} else if (genome_build == "hg19") {
  myc_exon_3_midpoint <- 128752921
  loci <- list(
    IGK = list(chrom = "chr2", start = 88299037, end = 90860939),
    MYC = list(chrom = "chr8", start = 127405427, end = 131774392),
    IGH = list(chrom = "chr14", start = 105055976, end = 106983247),
    IGL = list(chrom = "chr22", start = 21385754, end = 24247719))
} else {
  stop("Incorrect genome build")
}


# Load and tidy data --------------------------------------------------------------------------

bedpe <- fread(paths$input_bedpe)

if ("#CHROM_A" %in% colnames(bedpe)) {
  bedpe <- rename(bedpe, CHROM_A = "#CHROM_A")
}

colnames(bedpe) <- tolower(colnames(bedpe))

is_chr_prefix <- any(paste0("chr", c(1:22, "X")) %in% bedpe$chrom_a)

if (is_chr_prefix) {
  chroms <- paste0("chr", c(1:22, "X"))
} else {
  chroms <- as.character(c(1:22, "X"))
  loci <- map(loci, ~ map_at(.x, "chrom", ~ sub("chr", "", .)))
}

ig_chr_map <- map(loci, "chrom") %>% set_names(names(.), .)


# Create BED file for UCGC Genome Browser -----------------------------------------------------

colours <- list(BND = "34,138,96", DUP = "248,68,31", DEL = "68,90,152",
                INS = "204,65,139", INV = "97,175,23" )

intrachromosomal <-
  bedpe[filter == "PASS"] %>%
  .[type != "BND"] %>%
  .[, list(chrom = chrom_a, start = start_a, end = end_b,
           paste0(tumour_id, "_", type),
           1000, ".", start2 = start_a, end2 = end_b, colours[type])]

interchromosomal <-
  bedpe[filter == "PASS"] %>%
  .[type == "BND"] %>%
  .[, list(chrom = chrom_a, start = start_a, end = start_a,
           paste0(tumour_id, "_", type, "_", chrom_b, ":", start_b),
           1000, ".", start2 = start_a, end2 = start_a, colours[type])]

interchromosomal <-
  bedpe[filter == "PASS"] %>%
  .[type == "BND"] %>%
  .[, list(chrom = chrom_b, start = start_b, end = start_b,
           paste0(tumour_id, "_", type, "_", chrom_a, ":", start_a),
           1000, ".", start2 = start_b, end2 = start_b, colours[type])] %>%
  rbind(interchromosomal, .)

bed <- rbind(intrachromosomal, interchromosomal)

bed <- bed[chrom %in% chroms]


# Tabulate MYC translocations -----------------------------------------------------------------

myc <- bedpe[, .(chrom_a, start_a, end_a, strand_a,
                 chrom_b, start_b, end_b, strand_b,
                 tumour_id, 1000, name_a, name_b)]

myc <- myc[
  (chrom_a == loci$MYC$chrom & start_a > loci$MYC$start &
     end_a < loci$MYC$end & chrom_b %in% names(ig_chr_map)) |
    (chrom_b == loci$MYC$chrom & start_b > loci$MYC$start &
       end_b < loci$MYC$end & chrom_a %in% names(ig_chr_map))]

myc <-
  myc[chrom_a != loci$MYC$chrom] %>%
  setnames(
    c("chrom_a", "start_a", "end_a", "strand_a", "name_a",
      "chrom_b", "start_b", "end_b", "strand_b", "name_b"),
    c("chrom_b", "start_b", "end_b", "strand_b", "name_b",
      "chrom_a", "start_a", "end_a", "strand_a", "name_a")) %>%
  rbind(myc[chrom_a == loci$MYC$chrom], .)

myc <- myc[
  chrom_b == loci$IGK$chrom & start_b > loci$IGK$start & start_b < loci$IGK$end |
    chrom_b == loci$IGH$chrom & start_b > loci$IGH$start & start_b < loci$IGH$end |
    chrom_b == loci$IGL$chrom & start_b > loci$IGL$start & start_b < loci$IGL$end ]

myc <- myc[(start_a < myc_exon_3_midpoint & strand_a == "-") |
             (start_a > myc_exon_3_midpoint & strand_a == "+")]

myc_colours <-
  colours[1:3] %>%
  set_names(map_chr(loci[c("IGH", "IGK", "IGL")], "chrom"))

myc_bed <- rbind(
  myc[, .(chrom = chrom_a, start = start_a, end = end_a,
          paste0(tumour_id, "_", ig_chr_map[chrom_b], "_", chrom_b, ":", start_b),
          1000, ".",
          start2 = start_a, end2 = end_a, myc_colours[chrom_b],
          ".", ".", ".", name_a)],
  myc[, .(chrom = chrom_b, start = start_b, end = end_b,
          paste0(tumour_id, "_", ig_chr_map[chrom_b], "_", chrom_a, ":", start_a),
          1000, ".",
          start2 = start_b, end2 = end_b, myc_colours[chrom_b],
          ".", ".", ".", name_a)])


# Output --------------------------------------------------------------------------------------

fwrite(bedpe[filter == "PASS"], paths$output_bedpe, sep = "\t")

header <- 'track name="All SVs" description="Manta SVs (All)" itemRgb="On"'
fwrite(list(header), paths$output_bed_pass, quote = FALSE)
fwrite(bed, paths$output_bed_pass, sep = "\t", col.names = FALSE, append = TRUE)

myc_header <- 'track name="IG-MYC SVs" description="Manta SVs (IG-MYC)" itemRgb="On"'
fwrite(list(myc_header), paths$output_bed_myc, quote = FALSE)
fwrite(myc_bed, paths$output_bed_myc, sep = "\t", col.names = FALSE, append = TRUE)
