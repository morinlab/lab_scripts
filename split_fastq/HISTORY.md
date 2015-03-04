Release History
===============

1.1.0
-----

**Features and Enhancements**

- When one chunk is created, the output FASTQ is simply symlinked to the input FASTQ

**Bugfixes**

- When only one FASTQ file is created, now outputs a chunk name in interval file
- Chunk names include absolute paths instead of relative
- Now properly handles either one or two input FASTQ files (regarding output file names)


1.0.3
-----

**Bugfixes**

- When FASTQ read count is less than num_reads, no longer creates chunk_0


1.0.2
-----

**Features and Enhancements**

- Added num_buffer option to prevent memory bloat (enabled by default)
- Added default values for num_reads and num_buffer


1.0.1
-----

**Features and Enhancements**

- Using output prefix instead of output directory
- Output interval file contains suffixes instead of numbers


1.0.0
-----

Initial version.
