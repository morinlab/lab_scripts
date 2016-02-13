Create a local Makefile that contains the following. 

```makefile
# Export variable values to recursive Make calls
.EXPORT_ALL_VARIABLES:
NUM_THREADS := 15
BAM_TYPE := exome
BAM_DIR := /projects/fl_dart/FL_Okosun_Fitzgibbon_2015/bam_files
SYMLINK_DIR := /projects/rmorin/projects/nhl_meta_analysis/data/bam_originals/FL_Okosun_Fitzgibbon_2015
EGA_CLIENT := /projects/rmorin/software/centos-6/ega_demo_client/2.1.5/EgaDemoClient.jar
EGA_LOGIN := /home/bgrande/.ega_login
ENCRYPT_KEY := ihopeegaworks

include /path/to/lab_scripts/Makefile
```
