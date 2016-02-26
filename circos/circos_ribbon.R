################ make bed file for circos plot
#############################3

args <- commandArgs(trailingOnly = TRUE)
#usage: Rscript this_script.R tabulated_events.txt chromosome_details.txt bin_length
#chromosome_details file needs to use the circos naming style (hs prefixes) and the tabulated_events file should have chr prefixes

input_file = args[1]

chromosome_details = args[2]

bin_length = args[3]

### chromosome endpoints:

hs.karyotype.data <- read.table(chromosome_details, stringsAsFactors=FALSE, header=FALSE,sep='\t');
colnames(hs.karyotype.data) <- c('chr', 'Position');

##### add mitochondria to list of chromosomes 

chrEnds <- matrix(c(hs.karyotype.data$Position, 16571), nrow = 1);
colnames(chrEnds) <- c(hs.karyotype.data$chr, 'hsM');

######### translocation data

### recurrent translocations
# file with format:  chr1 12334326 TRA  chr2 11212344 5  ('FromChr' Position type 'ToChr' Position #recurrences)

trx.data.all <- read.table(input_file, header=FALSE,sep='\t');



########### circos ribbons need large regions otherwise show up as lines if given only one position
########### adjust the endpoints of the start and end of the ribbons, centered on the actual position, adjusted by number samples
########

adjustendpoints <- function(trx, endpts) {

	#colnames(trx) <- c('FromChr', 'ToChr', 'type', 'FromPos', 'ToPos', 'NumberSamples')
	#print(trx)
	trx$FromChr <- gsub('chr', 'hs', trx$FromChr)
	trx$ToChr <- gsub('chr', 'hs', trx$ToChr)

	trx$RibbonStart1 <- trx$FromPos - (trx$NumberSamples*bin_length)/2
	trx$RibbonStart2 <- trx$FromPos + (trx$NumberSamples*bin_length)/2 
	#print(trx)
	#q()
	for (i in (1:length(trx$FromChr))) {
		if (trx[i,'NumberSamples'] != 1 & trx[i,'RibbonStart2'] > endpts[1,trx[i,'FromChr']]) {
			trx[i,'RibbonStart1'] <- trx[i,'RibbonStart1'] + (endpts[1,trx[i,'FromChr']] - trx[i,'RibbonStart2'])
			trx[i,'RibbonStart2'] <- endpts[1,trx[i,'FromChr']]
		}
	}

	for (i in (1:length(trx$FromChr))) {
		if (trx[i,'NumberSamples'] != 1 & trx[i,'RibbonStart1'] < 0) {
			trx[i,'RibbonStart2'] <- -1*(trx[i,'RibbonStart1']) + trx[i,'RibbonStart2']
			trx[i,'RibbonStart1'] <- 0
		}
	}

	#print(trx)
	trx$RibbonEnd1 <- trx$ToPos - (trx$NumberSamples*bin_length)/2
	trx$RibbonEnd2 <- trx$ToPos + (trx$NumberSamples*bin_length)/2 
	#print(trx)
	#3q()
	for (i in (1:length(trx$ToChr))) {
		if (trx[i,'NumberSamples'] != 1 & trx[i,'RibbonEnd2'] > endpts[1,trx[i,'ToChr']]) {
			trx[i,'RibbonEnd1'] <- trx[i,'RibbonEnd1'] + (endpts[1,trx[i,'ToChr']] - trx[i,'RibbonEnd2'])
			trx[i,'RibbonEnd2'] <- endpts[1,trx[i,'ToChr']]
		}
	}

	for (i in (1:length(trx$ToChr))) {
		if (trx[i,'NumberSamples'] != 1 & trx[i,'RibbonEnd1'] < 0) {
			trx[i,'RibbonEnd1'] <- 0
			trx[i,'RibbonEnd2'] <- -1*(trx[i,'RibbonEnd1']) + trx[i,'RibbonEnd2']
		}
	}
	print(trx)
	return(trx)
}

########## format data to make bedfile required by circos

make_circos_table <- function(trxdata) {
    #print(trxdata)
    #q()
	trx.circos.data <- subset(trxdata, trxdata$NumberSamples != 1)
	#FromChr, RibbonStart1, RibbonEnd1, ToChr, RibbonStart2, RibbonEnd2  
	trx.circos.table <- trx.circos.data[,c('FromChr', 'RibbonStart1', 'RibbonStart2', 'ToChr', 'RibbonEnd1', 'RibbonEnd2', "NumberSamples"  )]
	trx.circos.table$NumberSamples <- gsub("([0-9]+)","value=\\1",trx.circos.table$NumberSamples)
	trx.circos.table$NumberSamples = paste(trx.circos.table$NumberSamples,",type=",trx.circos.data$type,sep="")
	trx.circos.table$RibbonStart2 <- format(trx.circos.table$RibbonStart2, scientific=FALSE)
	return(trx.circos.table)
}




colnames(trx.data.all) <- c('FromChr', 'ToChr', 'type', 'FromPos', 'ToPos', 'NumberSamples')
all.trx <- adjustendpoints(trx.data.all, chrEnds);
all.circos.table <- make_circos_table(all.trx);
wd =getwd()
write.table(all.circos.table, file="ctx_circos.bed", sep='\t', quote=FALSE,row.names=FALSE, col.names=FALSE);

q()

