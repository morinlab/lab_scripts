#!/bin/bash

# test if a file provided, otherwise read from stdin
[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"

(>&2 echo "Generating MutationMapper input from MAF...")

# filter out synonymous variants and 
# arrange columns in mutation mapper format
cat $input | tee >(egrep -h "^(#|Hugo_Symbol)") \
	| awk 'BEGIN {FS=OFS="\t"} $9 ~ /Splice_Site|Nonsense_Mutation|Frame_Shift_Del|Frame_Shift_Ins|Nonstop_Mutation|Translation_Start_Site|In_Frame_Ins|In_Frame_Del|Missense_Mutation/' \
	| tail -n+3 \
	| awk 'BEGIN {FS=OFS="\t"} {print $1,$16,$37,$9,$5,$6,$7,$11,$13}' \
	| sort \
	| uniq \
	| sed '1 iHugo_Symbol	Sample_ID	Protein_Change	Mutation_Type	Chromosome	Start_Position	End_Position	Reference_Allele	Variant_Allele'
