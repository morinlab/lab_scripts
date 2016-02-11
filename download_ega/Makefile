# Assumes that the datasets are listed in a file called datasets.txt

# If a recipe fails, delete the output file if it exists (most likely truncated)
.DELETE_ON_ERROR:

# Keep all intermediate files (I manually perform clean-up)
.SECONDARY:

.PHONY: all
all: task.complete

.PHONY: clean
clean: 
	-rm -rf ${BAM_DIR} ${SYMLINK_DIR}/*.bam targets.txt *.files.txt task.complete

task.complete: targets.txt
	mkdir -p ${BAM_DIR} ${SYMLINK_DIR}
	${MAKE} `<$<`
	touch $@

targets.txt: all.files.txt
	cut -f3 $< | perl -pe "s:(.*):${SYMLINK_DIR}/\1.bam:" > $@

all.files.txt: datasets.txt
	FILES=$(addsuffix .files.txt, $(shell cat $<)); \
	${MAKE} $$FILES; \
	cat $$FILES > $@

%.files.txt: 
	java -jar ${EGA_CLIENT} \
			-pf ${EGA_LOGIN} \
			--listdatasetfiles $* \
		| grep ".bam" \
		| awk 'BEGIN {FS="  "; OFS="\t";} {gsub(/\//,"_",$$2); sid=$$2; sid=gensub(/_EGAR[0-9]+_([^.]+).*/,"\\1", "g", sid); print $$4,$$2,sid; }' \
	> $@

${SYMLINK_DIR}/%.bam: ${BAM_DIR}/.%.decrypted
	ln -s `<$<` $@

${BAM_DIR}/.%.decrypted: ${BAM_DIR}/.%.downloaded
	cd ${BAM_DIR}; \
	java -jar ${EGA_CLIENT} \
		-pf ${EGA_LOGIN} \
		--decryptionkey ${ENCRYPT_KEY} \
		--decrypt `<$<`
	perl "s:(.*\.bam).*:\1:" $< > $@

${BAM_DIR}/.%.downloaded: ${BAM_DIR}/.%.requested all.files.txt
	cd ${BAM_DIR}; \
	java -jar ${EGA_CLIENT} \
		-pf ${EGA_LOGIN} \
		--downloadrequest `<$<` \
		--numthreads ${NUM_THREADS}
	ls `awk 'BEGIN {FS="\t"} $$3 == $* {$$2=gensub(/(.*\.bam).*/, "\\1", "g", $$2); print "${BAM_DIR}/" $$2 ".*"}' all.files.txt` > $@

${BAM_DIR}/.%.requested: all.files.txt
	EGAF=`awk 'BEGIN {FS="\t"} $$3 == "$*" {print $$1}' $<`; \
	REQ=`echo $$EGAF | perl -pe "s:^:req_:"`; \
	java -jar ${EGA_CLIENT} \
		-pf ${EGA_LOGIN} \
		--requestbyfileid $$EGAF \
		--reencryptionkey ${ENCRYPT_KEY} \
		-label $$REQ; \
	echo $$REQ > $@
