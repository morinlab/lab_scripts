# Assumes that vdb-config is already configured

# If a recipe fails, delete the output file if it exists (most likely truncated)
.DELETE_ON_ERROR:

# Keep all intermediate files (I manually perform clean-up)
.SECONDARY:

.PHONY: all
all: task.complete

.PHONY: clean
clean: sra_cache.txt
	-rm -rf targets.txt ${BAM_DIR}/*.bam `<$<` ${SYMLINK_DIR}/*.bam task.complete
	-rm $<

.PHONY: clean_cache
clean_cache: sra_cache.txt
	-rm -f `<$<`
	-rm $<

task.complete: targets.txt
	mkdir -p ${BAM_DIR} ${SYMLINK_DIR} ${TMP_DIR} ${NCBI_DIR}
	${MAKE} `<$<`
	touch $@
	
targets.txt: SraRunTable.txt
	grep -v BioSample_s $< \
		| cut -f ${SAMPLE_ID_COL} \
		| sort \
		| uniq \
		| perl -pe "s:(.*):${SYMLINK_DIR}/\1.${BAM_TYPE}.bam:" \
	> $@

sra_cache.txt: SraRunTable.txt
	grep -v BioSample_s $< \
		| cut -f ${SRR_ID_COL} \
		| perl -pe "s:(.*):${NCBI_DIR}/sra/\1.sra.cache:" \
	> $@

${SYMLINK_DIR}/%.${BAM_TYPE}.bam: SraRunTable.txt
	BAM=`awk -v sid="$*" 'BEGIN { FS="\t"; } $$${SAMPLE_ID_COL} == sid { print $$${SRR_ID_COL}; }' $< | head -1 | perl -pe "s:(.*):${BAM_DIR}/\1.${BAM_TYPE}.bam:"`; \
	${MAKE} $$BAM; \
	ln -s $$BAM $@

${BAM_DIR}/%.${BAM_TYPE}.bam:
	cd ${NCBI_DIR}; \
	${SAMDUMP} --disable-multithreading --header --unaligned $* \
		| ${SAMTOOLS} view -b - \
	> $@
	${SAMTOOLS} flagstat $@ > $@.flagstat 2> $@.flagstat.log
	-rm ${NCBI_DIR}/sra/$*.sra.cache  # Post-download clean-up to prevent cache build-up

