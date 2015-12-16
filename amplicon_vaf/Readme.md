amplicon_vaf
=================
```
amplicon_vaf.py --maf {variants.maf} --bam {samples.bam} --ref {reference.fa}
```

Required MAF Fields
====================

| Chromosome | Start_Position | Variant_Type | Reference_Allele | Tumor_Seq_Allele1 | Tumor_Seq_Allele2 |
| ---------- | -------------- | ------------ | ---------------- | ----------------- | ----------------- |
|          3 |      178952085 |          SNP |                A |                 A |                 G |

Fields Appended to MAF in Output
================================

| Reference_Counts | Alternate_Counts |      VAF |                     Sample_ID |
| ---------------- | ---------------- | -------- | ----------------------------- |
|             4702 |             3272 | 0.410334 | <basename of .bam input file> |
