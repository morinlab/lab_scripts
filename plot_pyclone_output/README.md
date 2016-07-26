Replot PyClone's loci-level CCF plots (`plots/loci/scatter.png`) with mutations in genes of interest labelled
- Usage details: `$ Rscript plot_ccfs.R -h`  
- Pass sample MAFs using ``--mafs`` to filter labelled mutations to nonsilent variants only

Example:
```
$ Rscript /software/lab_scripts/.../plot_ccfs.R \
    PT001 PT001_Tumour,PT001_Plasma pyclone_working_dir/PT001 lymphoma_genes.txt \
    --mafs mafs/PT001_Tumour.maf,mafs/PT001_Plasma.maf
```
