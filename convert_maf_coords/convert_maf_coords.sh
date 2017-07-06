#!/bin/bash

# ----------------------------------------------------------------------------------------------- #
#
#   convert_maf_coords.sh
#   ~~~~~~~~~~~~~~~~~~~~~
#
#   Usage:
#       convert_maf_coords.sh input.maf hg38ToHg19.over.chain.gz output.maf [crossmap|liftover]
#
#   The default mode uses CrossMap.py to convert between genomic coordinate systems, because it
#   tends to drop fewer variants than liftOver. In one test, CrossMap dropped 3.7% of the
#   898,384 input variants compared to liftOver's 7.0% drop rate.
#
#   The script assumes that there are no pipe (|) characters in your MAF file. If this is not
#   the case, you can specify your own MAF column separator using the MAFCOLSEP variable.
#
#   By default, the script expects CrossMap.py or liftOver to be in your PATH. If that's not the
#   case, you can specify the program path using the CROSSMAP or LIFTOVER environment variables.
#
# ----------------------------------------------------------------------------------------------- #


set -euf -o pipefail

# Parameters
SCRIPT="$(basename $0)"
INPUT_MAF="$1"
INPUT_BED="$(mktemp /tmp/${SCRIPT}.input.XXXXXX)"
CHAIN_FILE="$2"
OUTPUT_BED="$(mktemp /tmp/${SCRIPT}.output.XXXXXX)"
OUTPUT_MAF="$3"
OUTPUT_UNMAPPED="${OUTPUT_MAF%.*}.unmapped.bed"
MODE="${4:-crossmap}"
LOG_FILE="${OUTPUT_MAF%.*}.log.txt"
MAFCOLSEP="${MAFCOLSEP:-|}"

# For colour
RED='\033[0;31m'
NC='\033[0m'

DATE=$(date --rfc-3339 seconds | cut -c1-19)
echo "[${DATE}] INFO: Storing log in ${LOG_FILE}." | tee "${LOG_FILE}" >&2

DATE=$(date --rfc-3339 seconds | cut -c1-19)
echo "[${DATE}] INFO: Checking whether there are \"${MAFCOLSEP}\" characters in the input" \
	"MAF file." | tee "${LOG_FILE}" >&2
if fgrep -q \'${MAFCOLSEP}\' "${INPUT_MAF}"; then
	echo -e "${RED}[${DATE}] ERROR: ${INPUT_MAF} contains \"${MAFCOLSEP}\" characters, which" \
		"are used to separate MAF columns in the intermediate BED files. Please set the" \
		"MAFCOLSEP environment variable to change this.${NC}" | tee -a "${LOG_FILE}" >&2; exit 1;
fi

DATE=$(date --rfc-3339 seconds | cut -c1-19)
echo "[${DATE}] INFO: Creating temporary BED file (${INPUT_BED})" \
	"from input MAF file." | tee -a "${LOG_FILE}" >&2
awk 'BEGIN {FS=OFS="\t"} \
	{chrom = $5; start = $6; end = $7} \
	$10 == "SNP" {start = start - 1; end = end} \
	$10 == "DEL" {start = start - 1; end = end} \
	$10 == "INS" {start = start; end = end - 1} \
	$0 !~ /^(#|Hugo_Symbol)/ {row=gensub(/\t/, "'"${MAFCOLSEP}"'", "g", $0); print chrom, start, end, row}' \
	"${INPUT_MAF}" > "${INPUT_BED}"

if [[ "${MODE}" == "crossmap" ]]; then
	CROSSMAP="${CROSSMAP:-CrossMap.py}"
	command -v "${CROSSMAP}" >/dev/null 2>&1 || { \
		echo -e "${RED}[${DATE}] ERROR: ${CROSSMAP} is not in the PATH or the file doesn't" \
			"exist.${NC}" | tee -a "${LOG_FILE}" >&2; exit 1; }
	DATE=$(date --rfc-3339 seconds | cut -c1-19)
	echo "[${DATE}] INFO: Converting genomic coordinates using ${CROSSMAP}." \
		| tee -a "${LOG_FILE}" >&2
	"${CROSSMAP}" bed "${CHAIN_FILE}" "${INPUT_BED}" "${OUTPUT_BED}" 2>>"${LOG_FILE}"
	mv -f "${OUTPUT_BED}.unmap" "${OUTPUT_UNMAPPED}"
else
	DATE=$(date --rfc-3339 seconds | cut -c1-19)
	LIFTOVER="${LIFTOVER:-liftOver}"
	command -v "${LIFTOVER}" >/dev/null 2>&1 || { \
		echo -e "${RED}[${DATE}] ERROR: ${LIFTOVER} is not in the PATH or the file doesn't" \
			"exist.${NC}" | tee -a "${LOG_FILE}" >&2; exit 1; }
	echo "[${DATE}] INFO: Converting genomic coordinates using ${LIFTOVER}." \
		| tee -a "${LOG_FILE}" >&2
	"${LIFTOVER}" "${INPUT_BED}" "${CHAIN_FILE}" \
		"${OUTPUT_BED}" "${OUTPUT_UNMAPPED}" 2>>"${LOG_FILE}"
fi

NUM_UNMAPPED=$(grep -cv "^#" "${OUTPUT_UNMAPPED}")
if [[ ${NUM_UNMAPPED} > 0 ]]; then
	DATE=$(date --rfc-3339 seconds | cut -c1-19)
	echo -e "${RED}[${DATE}] WARNING: ${NUM_UNMAPPED} variants were not properly mapped to the" \
		"target genomic coordinate system.${NC}" | tee -a "${LOG_FILE}" >&2
fi

DATE=$(date --rfc-3339 seconds | cut -c1-19)
echo "[${DATE}] INFO: Converting output BED file (${OUTPUT_BED}) to MAF format." \
	| tee -a "${LOG_FILE}" >&2
grep "^Hugo_Symbol" "${INPUT_MAF}" > "${OUTPUT_MAF}"
awk 'BEGIN {FS=OFS="\t"} \
		{chrom = $1; start = $2; end = $3} \
		{row=gensub(/'"${MAFCOLSEP}"'/, "\t", "g", $4)} \
		$13 = "SNP" {start = start + 1; end = end} \
		$13 = "DEL" {start = start + 1; end = end} \
		$13 = "INS" {start = start; end = end + 1} \
		{print chrom, start, end, row}' \
		"${OUTPUT_BED}" \
	| awk 'BEGIN {FS=OFS="\t"} {$8 = $1; $9 = $2; $10 = $3; print $0}' \
	| cut -f4- >> "${OUTPUT_MAF}"

DATE=$(date --rfc-3339 seconds | cut -c1-19)
echo "[${DATE}] INFO: Cleaning up temporary BED files." | tee -a "${LOG_FILE}" >&2
rm -f "${INPUT_BED}" "${OUTPUT_BED}"

DATE=$(date --rfc-3339 seconds | cut -c1-19)
echo -e "[${DATE}] INFO: Done! You can find the following output files:
    MAF File:          ${OUTPUT_MAF}
    Unmapped Variants: ${OUTPUT_UNMAPPED}
    Log File:          ${LOG_FILE}" | tee -a "${LOG_FILE}" >&2
