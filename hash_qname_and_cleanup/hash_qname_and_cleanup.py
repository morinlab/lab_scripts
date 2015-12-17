import argparse
import os
import pysam
import glob
import cPickle
import gzip

def main():
    args = parse_args()
    og_bam = args.original_bam
    new_bams = args.new_bam
    new_fastqs = args.new_fastq
    read_set = args.read_set

    bam_files = []
    fastq_files = []

    original_reads = set()
    new_reads = set()

    if og_bam == None and read_set == None:
        print "--original_bam or --read_set parameter needs to be set."
        return
    elif not og_bam == None and not read_set == None:
        print "--original_bam and --read_set parameter both set. Only use one."
        return

    if new_bams == None and new_fastqs == None:
        print "--new_bam or --new_fastq parameter needs to be set."
        return
    elif not new_bams == None and not new_fastqs == None:
        print "--new_bam and --new_fastq parameter both set. Only use one."
        return

    if not og_bam == None and read_set == None:
        # original bam parameter set and read_set parameter not set.

        samfile = pysam.AlignmentFile(os.path.abspath(og_bam), 'rb')

        for read in samfile.fetch(until_eof=True):
            if read.is_read1:
                original_reads.add(read.query_name + '/1')
            elif read.is_read2:
                original_reads.add(read.query_name + '/2')
            else:
                original_reads.add(read.query_name)

        read_set_out = open('read_names.txt', 'wb')
        cPickle.dump(original_reads, read_set_out, -1)

    elif og_bam == None and not read_set == None:
        # read set parameter set and original bam parameter not set.
        read_set_in = open(read_set, 'rb')
        original_reads = cPickle.load(read_set_in)

    # Parse new bam or new fastqs into set
    if not new_bams == None and new_fastqs == None:
        if len(new_bams) == 1:
            if '*' in new_bams[0]:
                bam_files = glob.glob(new_bams[0])
            else:
                bam_files.append(new_bams[0])
        else:
            bam_files = new_bams

        for bam_file in bam_files:
            samfile = pysam.AlignmentFile(os.path.abspath(bam_file), 'rb')

            for read in samfile.fetch(until_eof=True):
                if read.is_read1:
                    new_reads.add(read.query_name + '/1')
                elif read.is_read2:
                    new_reads.add(read.query_name + '/2')
                else:
                    new_reads.add(read.query_name)

    elif new_bams == None and not new_fastqs == None:
        if len(new_fastqs) == 1:
            if '*' in new_fastqs[0]:
                fastq_files = glob.glob(new_fastqs[0])
            else:
                fastq_files.append(new_fastqs[0])
        else:
            fastq_files = new_fastqs

        for fastq_file in fastq_files:
            fastq = None
            
            if '.gz' in fastq_file:
                fastq = gzip.open(os.path.abspath(fastq_file))
            else:
                fastq = open(os.path.abspath(fastq_file))

            i = 0
            for line in fastq:
                if i % 4 == 0:
                    qname = line.rstrip()
                    new_reads.add(qname[1:])
                i += 1

    # Check if the two sets are the same.
    original_reads = original_reads - new_reads

    if not len(original_reads) == 0:
        print original_reads
        raise ValueError("Reads are missing in the new files.")
    else:
        print "SUCCESS!"

    return


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--original_bam', nargs='?',
                        help='BAM that will be realigned.')
    parser.add_argument('--new_bam', nargs='*',
                        help='BAM file(s) that will undergo read counting to \
                        check for read integrity.')
    parser.add_argument('--new_fastq', nargs='*',
                        help='FASTQ file(s) that will undergo read counting to \
                        check for read integrity.')
    parser.add_argument('--read_set', nargs='?',
                        help='File containing reads of the original BAM.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
