#! /usr/bin/env Rscript

#This program is intended to have a file feeded into it and extracts the useful information out of it in order to generate readable and
#reliable plots for viewing. As of now the only files tested are seqz files.
#The program takes the raw data of the seqz file and from it takes the columns corresponding to: Chromosome number, chromosome position for
#the chrommosome in question, and the adjusted depth ratio. 
#Data from Whole exome sequencing, single chromosome from a WGS hve been tested thus far on this script.


#install.packages("ggplot2")
#install.packages("dplyr")

library(ggplot2)
library(dplyr)
library(optparse)

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-i", "--ivar"), action="store", default=NULL, type='character',
              help="the input Seqz file"),
  make_option(c("-o", "--ovar"), action="store", default=NULL, type='character',
              help="output file"))
 
  parser = OptionParser(option_list=option_list)
  opt = parse_args(parser)
#Above is a parser (optparse) that is used for command line options in this script
if (! c("ivar") %in% names(opt) | ! c("ovar") %in% names(opt)) {
  print_help(parser)
  q()
}
   #A seqz file is a tab separated text file with column headers. The file has 14 columns. The first 3 columns are derived from the original 
   #'pileup file'
   #fast is a logical, if TRUE this file will be pre-parsed to count the number of rows; can speed up the file reading on some systems
   #nrows is the number of rows to read from the file. Default is -1, which is all rows
   #takes a seqz files and arra nges the columns with their corresponding labels as either being characters or strings
read.seqz <- function (file, nrows = -1, fast = FALSE, gz = TRUE, header = TRUE,
    colClasses = c('character', 'integer', 'character', 'integer',
      'integer', 'numeric', 'numeric', 'numeric', 'character',
      'numeric', 'numeric', "character", "character", "character"),
    chr.name = NULL, n.lines = NULL, ...) {
	#will read a seqz or acgt file into R
   if (!is.null(n.lines) & is.null(chr.name)) fast <-  FALSE 
   if(fast && nrows == -1) {
	#fast is a logical, if TRUE this file will be pre-parsed to count the number of rows; can speed up the file reading on some systems
	#nrows is the number of rows to read from the file. Default is -1, which is all rows
	#takes a seqz files and arra nges the columns with their corresponding labels as either being characters or strings
    if(gz) {
    	#gz is a logical, if TRUE will expect a gzipped file (the default)
       if (!is.null(chr.name)) {
          wc <- system(paste('gzip -d -c ',file,' | grep -c "^', chr.name, '\t"', sep = ''), intern = TRUE)
       } else {
          wc <- system(paste('gzip -d -c', file, '| wc'), intern = TRUE)
       }
    } else {
       if (!is.null(chr.name)) {
       	#if chr.name is specified, only the selected chromosome will be extracted instead of the entire file. 
          wc <- system(paste(paste('grep -c "^', chr.name, '\t"', sep = ''), file, sep = ' '), intern = TRUE)
       } else {
          wc <- system(paste('wc', file), intern = TRUE)
       }
    }
    if (is.null(chr.name)) {
       wc <- sub("^ +", "", wc)
       wc <- strsplit(wc, ' ')[[1]][1]
    }
    nrows <- max(as.integer(wc), 1)
    message('Reading ', nrows, ' lines...')
  }
   if (!is.null(chr.name)) {
      if (gz) {
         grep.part <- paste("gzip -d -c ", file," | grep '^", chr.name, "\t'", sep = "")
      } else {
         grep.part <- paste("grep '^", chr.name, "\t'", sep = "")
      }
      seqz.data   <- read.delim(pipe(grep.part), nrows = nrows, colClasses = colClasses, header = FALSE, ...)
      if (header == TRUE) {
      	#header is a logical indicating whether the file contains the names of the variables as its first line
         head       <- colnames(read.table(file, header = TRUE, nrows = 1 ))
         colnames(seqz.data) <- head
      }
      seqz.data
   } else {
      if (!is.null(n.lines)){
      	#vector of length 2 specifying the first and last line to read from the file. If specified, only the selected portion of the file
      	#will be used. Requires the "sed" UNIX utility.
         if (!is.numeric(n.lines) | length(n.lines) != 2) stop("n.lines must be a vector of 2 integers")
         n.lines <- round(sort(n.lines), 0)
         if (header == TRUE) {
            n.lines <- n.lines + 1
         }
         if(gz) {
            seqz.data <- read.delim(pipe(paste("gzip -d -c", file,"| sed -n '", paste(n.lines[1], n.lines[2], sep = ","),"p'")),
                       colClasses = colClasses, nrows = 1 + n.lines[2] - n.lines[1], header = FALSE,...)
         }  else{
            seqz.data <- read.delim(pipe(paste("sed -n '", paste(n.lines[1], n.lines[2], sep = ","),"p'", file)),
                       colClasses = colClasses, nrows = 1 + n.lines[2] - n.lines[1], header = FALSE, ...)
            #colClasses represents a character. A vector of classes to be assumed for the columns. By default the acgt and seqz file is expected
         }
         if (header == TRUE) {
            head  <- colnames(read.table(file, header = TRUE, nrows = 1 ))
            colnames(seqz.data) <- head
         }
         seqz.data
      } else {
         read.delim(file, nrows = nrows, colClasses = colClasses, header = header, ...)
      }
   }
}

gc.norm <- function (x, gc) {
   #gc.norm detects the bias in 'x' driven by 'gc'. Specifically for each value of 'gc', summary statistics are calculated
   #for the corresponding values of 'x'. These statistics can then be used to normalize 'x' for 'gc'.
   #x is the vector of values to be normalized by GC-content, typically depth ratio values.
   #gc.norm function is performed in the columns depth ratio and gc.percent, which were taken from the seqz file.
   dr.by.gc <- split(x = x, f = gc)
   raw <- t(sapply(dr.by.gc, quantile, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
   dr.by.gc.mean <- sapply(dr.by.gc, FUN = mean, na.rm = TRUE)
   dr.by.gc.median <- sapply(dr.by.gc, FUN = median, na.rm = TRUE)
   adj <- sweep(raw, 1, dr.by.gc.median, '/')
   list(raw = raw, adj = adj, gc.values = as.numeric(names(dr.by.gc)),
        raw.mean = dr.by.gc.mean, raw.median = dr.by.gc.median)
}
#gc.norm takes the gc values from the seqz file and does thr normalization of them. The raw gc values are divided by the median values which
#yeild the three quartiles of 0.25, 0.5 and 0.75. The 0.75 is normalized to a value of 1 as it is the median devided by itself.

gc.sample.stats <- function (file, gz = TRUE) { #if true (default) the function expects a gzipped file
	#gc.sample.stats is a function that will take the gzipped seqz file and extract from it and paste together the desired columns.
	#functions gc.norm and gc.stats are both contained in withing gc.sample.stats.
	#gc.sample.stats detects the bias in the depth ratio values driven by varying gc content.
   colClasses = c('character', 'numeric', 'numeric')
   if (gz) {
      seqz.data <- read.delim(pipe(paste('gzip -d -c', file, '| cut -f 1,6,10')), colClasses = colClasses)
      #unix pipe that takes columns 1,6,and,10 from the seqz file and joins them together (can enter any column number desired)
   } else {
      seqz.data <- read.delim(pipe(paste('cut -f 1,6,10', file)), colClasses = colClasses)
   }
   gc.stats <- gc.norm(x  = seqz.data$depth.ratio,
                       gc = seqz.data$GC.percent)
   #gc.norm detects the bias in 'x' driven by 'gc'. Specifically for each value of 'gc', summary statistics are calculated
   #for the corresponding values of 'x'. These statistics can then be used to normalize 'x' for 'gc'.
   #x is the vector of values to be normalized by GC-content, typically depth ratio values.
   #gc.norm function is performed in the columns depth ratio and gc.percent, which were taken from the seqz file.
   chr.ord  <- unique(seqz.data$chromosome)
   #the unique funtion will return an object if the same class as x, which is a data frame like object. Chromsome column is taken from seqz
   chr.dim  <- lapply(X = split(seqz.data$chromosome, seqz.data$chromosome), FUN = length)
   #lapply returns list same length as 'X', each element has FUN applied to it
   chr.dim  <- data.frame(chr = chr.ord, n.lines = do.call(rbind,chr.dim[chr.ord]))
   chr.dim$start <- cumsum(c(1, chr.dim$n.lines[-length(chr.dim$n.lines)]))
   chr.dim$end   <- chr.dim$start + chr.dim$n.lines - 1
   #Lists the start and end coordinates of the chromosomes.
   gc.stats$file.metrics <- chr.dim
   gc.stats
   #the output of the function gc.sample.stats is saved and assigned to gc.stats
}


gc.compensate <- function(seqzFile, output_file, normalization.method = "median") {
  #gc.compensate will take whole exome sequencing data and from the respective seqz file fed into it retreieve the columns corresponding to
  #the specific genomic coordinates as well as the compensated depth ratio after the gc correction is applied to the data.

  gc.stats = gc.sample.stats(seqzFile)
  chr.vect <- as.character(gc.stats$file.metrics$chr)
   if (normalization.method != "mean") {
      gc.vect  <- setNames(gc.stats$raw.median, gc.stats$gc.values)
   } else {
      gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
   }
  seqz.data = read.seqz(seqzFile, gz=TRUE) 
  # print(str(gc.compensate))
  seqz.data$adjusted.ratio = round(seqz.data$depth.ratio/gc.vect[as.character(seqz.data$GC.percent)], 3)
  write.table(seqz.data, sep="\t", file="TestSeqzNorm.tsv", row.names=FALSE)
  graph.ggplot2(seqz.data, output_file)
}

graph.ggplot2 <- function(seqz.data, output_file) {
 seqz.data %>% sample_frac(0.001) %>% ggplot(aes(x = position, y = adjusted.ratio, color = chromosome)) + geom_point() + facet_wrap(~chromosome)
 ggsave(output_file, plot = last_plot(), device = "pdf", dpi = 300, limitsize = TRUE)
#graph.ggplot2 is a very simple function that will take the outputted file from gc.compensate, after the cooumns have been appropriately selected and the 
#compensation correction is performed in the adjusted depth ratio and then proceed to plot the data points via ggplot2.
#dplyr should be installed as indictaed at the top level of the code in order to make use of the functions withing the package. 
#seqz.data will be piped to a function that will subset the data and then carry on to graph the genomic coordinate against the adjusted depth
}


seqzFile = opt$i 
output_file = opt$o
gc.compensate(seqzFile, output_file)





