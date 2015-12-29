import argparse
import pysam
import os
import glob
import gzip

"""
Description: Checks the integrity of FASTQ or BAM file(s)
compared to an original BAM file. Checking entails finding
the query name and position in pair of a read and applying
a hash function to produce an integer for the read. These
integers are summed producing a unique hash sum for the BAM
or FASTQ file(s). These hash sums can be compared to
determine if there is a difference between two sequence
files. This script also identifies supplementary alignments,
unpaired reads and split reads and ignores them when
calculating the hash sum.

Parameters:
--original_bam and --hash_sum_outfile parameters must be set together.
--new_fastqs and --hash_sum_infile or --new_bams and --hash_sum_infile
must be set together.
"""

def main():
    args = parse_args()
    original_bam = args.original_bam
    hash_sum_outfile = args.hash_sum_outfile
    hash_sum_infile = args.hash_sum_infile
    new_fastqs = args.new_fastqs
    new_bams = args.new_bams

    if original_bam and hash_sum_outfile:
        original_hash_sum, paired_reads = sum_bam(original_bam)
        write_hash_sum(original_hash_sum, hash_sum_outfile)

    elif hash_sum_infile and ( new_fastqs or new_bams ):
        old_hash_sum = long(hash_sum_infile.readline().rstrip())
        hash_sum_infile.close()

        if new_fastqs:
            new_hash_sum = sum_new_fastqs(new_fastqs)
        elif new_bams:
            new_hash_sum = sum_new_bams(new_bams)

        if not old_hash_sum == new_hash_sum:
            print new_hash_sum
            raise ValueError('New hash sum does not match original hash sum.')
    else:
        raise ValueError('Parameter error.')

    return

def write_hash_sum(hash_sum, hash_sum_outfile):
    hash_sum_outfile.write(str(hash_sum) + '\n')
    hash_sum_outfile.close()
    return

def sum_new_bams(new_bams):
    bams = None
    hash_sum = 0
    paired_reads = {}

    if len(new_bams) == 1:
        if '*' in new_bams[0]:
            bams = glob.glob(os.path.abspath(new_bams[0]))
        else:
            bams = [ new_bams[0] ]
    else:
        bams = [ os.path.abspath(b) for b in new_bams ]

    for bam in bams:
        hash_sum, paired_reads = sum_bam(bam, hash_sum, paired_reads)

    return hash_sum

def sum_bam(bam, hash_sum=0, paired_reads={}):
    sam = pysam.AlignmentFile(bam, 'rb')

    reads = sam.fetch(until_eof=True)

    for read in reads:

        if read.is_secondary:
            continue

        if read.is_supplementary:
            continue

        qname = read.query_name

        if qname not in paired_reads:
            paired_reads[qname] = [False, False]

        if read.is_read1:
            paired_reads[qname][0] = True

        elif read.is_read2:
            paired_reads[qname][1] = True

        else:
            hash_sum += hash(qname)
            continue
 
        if all(paired_reads[qname]):
            hash_sum += hash(qname + '/1')
            hash_sum += hash(qname + '/2')
            del paired_reads[qname]

    for qname in paired_reads.keys():
        hash_sum += hash(qname)

    sam.close()

    return hash_sum, paired_reads

def sum_new_fastqs(new_fastqs):
    fastqs = None
    hash_sum = 0

    if len(new_fastqs) == 1:
        if '*' in new_fastqs[0]:
            fastqs = glob.glob(os.path.abspath(new_fastqs[0]))
        else:
            fastqs = [ new_fastqs[0] ]
    else:
        fastqs = [ os.path.abspath(f) for f in new_fastqs ]

    for fastq in fastqs:
        fq = None

        if '.gz' in fastq:
            fq = gzip.open(fastq, 'rb')
        else:
            fq = open(fastq, 'r')

        i = 0
        for line in fq:
            if i % 4 == 0:
                qname = line.rstrip()[1:]
                hash_sum += hash(qname)
            i += 1

        fq.close()

    return hash_sum

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--original_bam',
                        help='Path to original BAM.')
    parser.add_argument('--hash_sum_outfile', type=argparse.FileType('w'),
                        help='Output file for the hash sum of original BAM reads.')
    parser.add_argument('--hash_sum_infile', type=argparse.FileType('r'),
                        help='File containing a hash sum that will be used to check \
                        integrity of either new fastq or new bam files.')
    parser.add_argument('--new_fastqs', nargs='*',
                        help='Path(s) or a single glob to FASTQ files to check for \
                        read name integrity.')
    parser.add_argument('--new_bams', nargs='*',
                        help='Path(s) or a single glob to FASTQ files to check for \
                        read name integrity.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()
