Morin Lab Scripts
=================

This repository serves to consolidate scripts for shared use in the Morin lab. Each script is housed in a separate subdirectory, which contains all related files, if applicable. 

To reduce redundant code, commonly used methods are added to a module in the `modules` subdirectory. Accordingly, users should add the path of this `modules` directory to their `PYTHONPATH` environment variable. 

    export PYTHONPATH=/path/to/lab_scripts/modules:$PYTHONPATH

A quick summary and example command for each script are included below. 

calc_vaf_db_strand
------------------

### Description

calc_vaf_db_strand.py calculates the variant allele fraction (VAF) at a given position while only considering read pairs for which the forward and reverse strand agree on the base call at the given position. Effectively, this suppresses errors due to incorrect base calling. This script was originally meant to calculate the VAF in amplicon sequencing data derived from circulating tumour DNA (ctDNA), where the forward and reverse strands are designed to overlap. 

    python calc_vaf_db_strand.py --in_bam input.sorted.bam --out_bam output.sorted.db_strand.bam --chromosome chr3 --position 38182641 --ref_allele T --mut_allele C

calc_fusion_vaf
---------------

### Requirements

- Pysam 0.8.0 (for now, pysam 0.8.1 does not work)

### Description

calc_fusion_vaf.py calculates the variant allele fraction (VAF) of fusions called by Factera. It does so by appending the fusion sequences to the reference genome and realigning the reads to this new reference. Then, reads supporting the fusion and wild-type alleles are counted and a VAF is calculated. Moreover, if a certain gene was targeted like in a capture-based approach, there's the option of only considering fusions involving this gene (_e.g._ `--gene EWSR1`).

    python calc_fusion_vaf.py --output_dir ./fusion_vaf/ --threads 4 hg19.fa ewings_sarcoma.R1.fastq.gz ewings_sarcoma.R2.fastq.gz ./factera_output/ewings_sarcoma.factera.fusions.txt
