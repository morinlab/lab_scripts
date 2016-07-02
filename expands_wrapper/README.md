
# EXPANDS wrapper

## Usage

```
$ Rscript /path/to/run_analysis.R seg_file seg_input_mode maf_file sample_name output_dir \  
  [--loh {0, 1, 2, 3}] \  
  [--max_score MAX_SCORE] \  
  [--precision PRECISION] \  
  [--cn_style {1, 2}] \  
  [--pyclone_dir PYCLONE_DIR] \  
  [--pyclone_only {TRUE, FALSE}]
```


```
positional arguments:
  seg                   Input segments file
  input_mode			Type of seg file: S (Sequenza), I (IGV-friendly seg file), T (Titan), O (augmented OncoSNP file)
  maf			        Input MAF file
  sample		    	Sample ID
  output_dir			Path to directory where all output will be saved

flags:
  -h, --help			show this help message and exit

optional arguments:
  -x, --opts OPTS		          RDS file containing argument values
  -l, --loh LOH		  	          0: ignore LOH events,
                                  1: include all copy-neutral LOH segments and their BAF in clustering; can help resolve clonal       clusters with few mutations, (recommended)
                                  2: include deletion LOH only,
                                  3: include all LOH [default: 1]
  -m, --max_score MAX_SCORE		max_score for EXPANDS [default: 2.25]
  -p, --precision PRECISION		precision for EXPANDS [default: 0.05]
  -c, --cn_style CN_STYLE			  1 for integer values,
                                      2 for rational numbers calculated from CN log ratios (recommended) [default: 2]
  --pyclone_dir PYCLONE_DIR			Specific output directory for PyClone files
  --pyclone_only PYCLONE_ONLY		TRUE: Generate PyClone input only, skip EXPANDS [default: FALSE]
```
