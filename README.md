Morin Lab Scripts
=================

This repository serves to consolidate scripts for shared use in the Morin lab. Each script is housed in a separate subdirectory, which contains all related files, if applicable. 

To reduce redundant code, commonly used methods are added to a module in the `modules` subdirectory. Accordingly, users should add the path of this `modules` directory to their `PYTHONPATH` environment variable. 

    export PYTHONPATH=/path/to/lab_scripts/modules:$PYTHONPATH

A quick summary and example command for each script are included below. 

calc_vaf_db_strand
------------------

calc_vaf_db_strand.py calculates the variant allele fraction (VAF) at a given position while only considering read pairs for which the forward and reverse strand agree on the base call at the given position. Effectively, this suppresses errors due to incorrect base calling. 

    python calc_vaf_db_strand.py --in_bam input.sorted.bam --out_bam output.sorted.db_strand.bam --chromosome chr3 --position 38182641 --ref_allele T --mut_allele C
