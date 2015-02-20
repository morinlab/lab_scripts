# calc_fold_enrichment

This script calculates useful statistics for capture-based experiments. As input, it takes a BAM alignment file and a set of intervals (BED file). It then produces a simple TSV file describing various statistics. Currently, calc_capture_stats calculates the following statistics:

- Average genome coverage
- Average target coverage
- Percent on-target reads
- Percent fold enrichment

### Requirements

- pysam (v0.8.1 or later)
- [cancer_api|https://github.com/brunogrande/cancer_api]
