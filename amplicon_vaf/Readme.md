amplicon_vaf
=================
```
amplicon_vaf.py --maf {variants.maf} --bam {samples.bam} --ref {reference.fa}
```

Required MAF Fields
====================
Below is an example .maf input file. The header line is required. Additional fields may be included, and the order of fields may be altered. Refer to the [MAF Specification](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification) for more details on MAF format.

| Chromosome | Start_Position | Variant_Type | Reference_Allele | Tumor_Seq_Allele1 | Tumor_Seq_Allele2 |
| ---------- | -------------- | ------------ | ---------------- | ----------------- | ----------------- |
|          3 |      178952085 |          SNP |                A |                 A |                 G |
|          5 |      112173673 |          DEL |               TC |                TC |                 - |
|          5 |      112176020 |          INS |                - |                 - |                AT |

An optional `Variant_ID` field may be included to help with downstream analysis with maf2matrix.R.

Fields Appended to MAF in Output
================================

| Reference_Counts | Alternate_Counts |      VAF |                     Sample_ID |
| ---------------- | ---------------- | -------- | ----------------------------- |
|             4702 |             3272 | 0.410334 | (basename of .bam input file) |

Running the maf2matrix.R Helper Script
======================================
The `maf2matrix.R` script will convert the output of `amplicon_vaf.py` into a matrix (variants x samples) of VAF values.
If the optional `Variant_ID` field is included, then that will be used for the row label.  If no `Variant_ID` is provided
then IDs will be constructed based on position and base change (eg: `3:178952085A>G`).
```
Rscript maf2matrix.R {amplicon_vaf_output.maf}
```
