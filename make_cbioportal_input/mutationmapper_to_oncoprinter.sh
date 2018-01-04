	#!/bin/bash

# test if a file provided, otherwise read from stdin
[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"

(>&2 echo "Generating OncoPrinter input from MutationMapper input...")

# assumes a list of sample IDs exists (samples.txt)
# convert mutation mapper format to input format
# for oncoprinter
cat $input | tail -n+2 \
	| awk 'BEGIN {FS=OFS="\t"} {$3 = gensub(/p\./, "", "g", $3)} \
		$3 ~ /splice|fs|\*|\?/ {$4 = "TRUNC"} $3 ~ /del|>/ {$4 = "INFRAME"} \
		$3 !~ /splice|fs|\*|del|>/ {$4 = "MISSENSE"} {print $2,$1,$3,$4}' \
	| sort -k2,2 -k1,1 \
	| uniq \
	| sed '1 iSample	Gene	Alteration	Type'
