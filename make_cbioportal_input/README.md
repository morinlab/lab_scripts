
# Usage

With a set of MAF files:

```
$ ./concat_mafs.sh *.maf | ./maf_to_mutationmapper.sh | tee mutationmapper.txt \
| ./mutationmapper_to_oncoprinter.sh > oncoprinter.txt
```


With a single concatenated MAF file:

```
$ ./maf_to_mutationmapper.sh all.maf | tee mutationmapper.txt | ./mutationmapper_to_oncoprinter.sh > oncoprinter.txt
```