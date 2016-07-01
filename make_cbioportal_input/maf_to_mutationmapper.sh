#!/bin/bash

# test if a file provided, otherwise read from stdin
[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"

# filter out synonymous variants and 
# arrange columns in mutation mapper format
cat $input | tee >(egrep -h "^(#|Hugo_Symbol)") \
	| awk 'BEGIN {FS=OFS="\t"} $37 ~ /p\.[^=]/ { print }' \
	| awk 'BEGIN {FS=OFS="\t"} $51 !~ /synonymous_variant/ { print }' \
	| tail -n+3 \
	| awk 'BEGIN {FS=OFS="\t"} {print $1,$16,$37,$9,$5,$6,$7,$11,$13}' \
	| sort \
	| uniq \
	| sed '1 iHugo_Symbol	Sample_ID	Protein_Change	Mutation_Type	Chromosome	Start_Position	End_Position	Reference_Allele	Variant_Allele'