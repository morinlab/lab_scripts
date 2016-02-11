Create a local Makefile that contains the following. 

```makefile
# Export variable values to recursive Make calls
.EXPORT_ALL_VARIABLES:
SAMPLE_ID_COL := 5
SRR_ID_COL := 11
BAM_TYPE := exome
BAM_DIR := /projects/lymphoma_dart/FL_DLBCL_Pasqualucci_DallaFavera_2014/bam_files
SYMLINK_DIR := /projects/rmorin/projects/nhl_meta_analysis/data/bam_originals/FL_DLBCL_Pasqualucci_DallaFavera_2014
NCBI_DIR := /projects/bgrande/ncbi
SAMTOOLS := /projects/rmorin/software/centos-6/samtools/1.2/bin/samtools
SAMDUMP := /projects/rmorin/software/centos-6/sratoolkit/2.5.1/bin/sam-dump

include /path/to/lab_scripts/Makefile
```
