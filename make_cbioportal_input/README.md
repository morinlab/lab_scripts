
# Usage

With a set of MAF files:

```
$ ./concat_mafs.sh *.maf | tee all.maf \
| ./maf_to_mutationmapper.sh | tee mutationmapper.txt \
| ./mutationmapper_to_oncoprinter.sh > oncoprinter.txt
```


With a single concatenated MAF file:

```
$ ./maf_to_mutationmapper.sh all.maf | tee mutationmapper.txt \
| ./mutationmapper_to_oncoprinter.sh > oncoprinter.txt
```

Note:
If there are samples in your cohort with no mutations, it's a good idea to append the list of sample IDs to the OncoPrinter input (e.g. `cat samples.txt >> oncoprinter.txt`) so that the OncoPrinter alteration percentages are calculated correctly.