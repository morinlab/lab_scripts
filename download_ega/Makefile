# Assumes that the datasets are listed in a file called datasets.txt

# If a recipe fails, delete the output file if it exists (most likely truncated)
.DELETE_ON_ERROR:

# Keep all intermediate files (I manually perform clean-up)
.SECONDARY:

# Recipe for download, checking and symlinking all of the BAM files
.PHONY: all
all: task.complete

# Recipe for deleting output files
.PHONY: clean
clean: 
	-rm -rf ${BAM_DIR} ${SYMLINK_DIR}/*.bam targets.txt *.files.txt task.complete

# Ensure that the directories exist and make the target symlinks
task.complete: targets.txt
	mkdir -p ${BAM_DIR} ${SYMLINK_DIR}
	${MAKE} `<$<`
	touch $@

# Extract the sample name and and create the target symlink
.PHONY: targets.txt
targets.txt: all.files.txt
	cut -f3 $< | perl -pe "s:(.*):${SYMLINK_DIR}/\1.${BAM_TYPE}.bam:" > $@

# Create a list of BAM files for each dataset and merge them together
all.files.txt: datasets.txt
	FILES=$(addsuffix .files.txt, $(shell cat $<)); \
	${MAKE} $$FILES; \
	cat $$FILES > $@

# Create a list of BAM files for a given dataset
%.files.txt: 
	java -jar ${EGA_CLIENT} \
			-pf ${EGA_LOGIN} \
			--listdatasetfiles $* \
		| grep ".bam" \
		| grep -v ".bai" \
		| awk 'BEGIN {FS="  "; OFS="\t";} {gsub(/\//,"_",$$2); sid=$$2; sid=gensub(/_EGAR[0-9]+_([^.]+).*/,"\\1", "g", sid); print $$4,$$2,sid; }' \
	> $@

# Symlink the BAM to the symlink directory
${SYMLINK_DIR}/%.${BAM_TYPE}.bam: ${BAM_DIR}/.%.decrypted ${BAM_DIR}/.%.flagstat
	ln -fs `<$<` $@

# Run flagstat on the BAM file to check its integrity
${BAM_DIR}/.%.flagstat: ${BAM_DIR}/.%.decrypted
	${SAMTOOLS} flagstat `<$<` > `<$<`.flagstat 2> `<$<`.flagstat.log
	-rm -f `<$<`.flagstat.log  # Only remove the log if the previous command worked
	touch $@

# Decrypt the downloaded BAM file and store the name of the BAM file in the target file
${BAM_DIR}/.%.decrypted: ${BAM_DIR}/.%.downloaded
	cd ${BAM_DIR}; \
	java -jar ${EGA_CLIENT} \
		-pf ${EGA_LOGIN} \
		--decryptionkey ${ENCRYPT_KEY} \
		--decrypt `<$<`
	perl -pe "s:(.*\.bam).*:\1:" $< > $@

# Download the encrypted BAM file and store the list of encrypted files in the target file
${BAM_DIR}/.%.downloaded: ${BAM_DIR}/.%.requested 
	cd ${BAM_DIR}; \
	java -jar ${EGA_CLIENT} \
		-pf ${EGA_LOGIN} \
		--downloadrequest `<$<` \
		--numthreads ${NUM_THREADS}
	ls $(shell awk 'BEGIN {FS="\t"} $$3 == "$*" {$$2=gensub(/(.*\.bam).*/, "\\1", "g", $$2); print "${BAM_DIR}/" $$2 ".*"}' all.files.txt) > $@

# Make a request for a file and store the request ID in the target file
${BAM_DIR}/.%.requested: 
	EGAF=`awk 'BEGIN {FS="\t"} $$3 == "$*" {print $$1}' all.files.txt`; \
	REQ=`echo $$EGAF | perl -pe "s:^:req_:"`; \
	java -jar ${EGA_CLIENT} \
		-pf ${EGA_LOGIN} \
		--requestbyfileid $$EGAF \
		--reencryptionkey ${ENCRYPT_KEY} \
		-label $$REQ; \
	echo $$REQ > $@
