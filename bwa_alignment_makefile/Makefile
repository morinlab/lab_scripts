# Thanks to Biostars.org user 'Ketil' (https://www.biostars.org/u/1323/)
# https://www.biostars.org/p/9472/#9486

##################################################
# Align illumina reads

%.bam: %_L001_R1_001.fastq.gz %_L001_R2_001.fastq.gz 
	bwa mem -t 4 $(GENOME) $*_L001_R1_001.fastq.gz $*_L001_R2_001.fastq.gz | samtools view -bS -o $@ -

##################################################
# Convert to sorted and indexed BAM file

%.sorted.bam: %.bam
	samtools sort -@ 4 $< $*.sorted
	rm $<

%.sorted.bam.bai: %.sorted.bam
	samtools index $*.sorted.bam

%.sorted.rmdup.bam: %.sorted.bam
	samtools rmdup $*.sorted.bam $*.sorted.rmdup.bam

%.sorted.rmdup.bam.bai: %.sorted.rmdup.bam
	samtools index $*.sorted.rmdup.bam

.PRECIOUS: %.sorted.bam %.sorted.rmdup.bam
