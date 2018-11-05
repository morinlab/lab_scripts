# get_bams.py

This script uses the GSC's internal REST API to retrieve the BAM file paths for given library IDs. 
For more information, read the [script](get_bams.py) header. You may also get the command-line interface by running: 

```console
$ python ./get_bams --help
```


## Dependencies

This script has no software dependencies other than Python 3 and the requests HTTP package. 
You will need to create an INI configuration file containing your GIN credentials. 
See the [script](get_bams.py) header for more details.


## Usage

This script is simple to use. 
Below is an example command, where I'm retrieving the paths for merged genome BAM files using the GSC API. 
The file paths are outputted to stdout by default in TSV format.

Here, you can notice a warning on stderr that multiple BAM files were returned for at least one library ID. 
It turns out that these libraries were aligned to GRCh37 and GRCh38, causing the warning. 
By default, get_bams.py will return the latest BAM file using the timestamps returned by the API.

```console
$ python ./get_bams -t genome A47817 A47818 A57139
A47817  /path/to/A47817.bam
A47818  /path/to/A47818.bam
A57139  /path/to/A57139.bam
WARNING: At least one library ID returned multiple BAM files. The most recent one was returned by default. To obtain all of the BAM file paths, enable -a. To inspect the API results more closely, enable -v or -vv. To disable this warning, enable -q. For more information about these options, run with --help.
```
