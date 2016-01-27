import argparse
import pysam
import logging
import subprocess
import glob
import os

def main():
    args = parse_args()
    bam = args.bam
    output_dir = args.output_dir
    interval_file = args.interval_file
    num_reads = args.num_reads

    #if not os.path.exists(output_dir):
    #    os.path.mkdirs(output_dir, 0775)

    logging.basicConfig(filename='split_reads.log', level=logging.INFO,
                        filemode='w')

    sam = pysam.AlignmentFile(bam, 'rb')

    reads = sam.fetch(until_eof=True)

    paired_dict = {}

    chunk = 0

    paired_string = 'paired_{}.fastq.gz'.format(chunk)
    paired_process = spawn_gzip(os.path.join(output_dir, paired_string)) 

    chunk += 1
    
    unpaired_string = 'unpaired_{}.fastq.gz'.format(chunk)
    unpaired_process = spawn_gzip(os.path.join(output_dir, unpaired_string))

    read_count = 0
    unpaired_read_count = 0

    chunk_files = [ paired_string, unpaired_string ]

    for read in reads:

        if read.is_secondary:
            continue

        #if read.has_tag('SA') and not read.is_supplementary:
        #    logging.info(str(read))

        if read.is_supplementary:
            continue

        if read.is_duplicate:
            continue

        if read.is_qcfail:
            continue

        qname = read.query_name

        seq = read.query_sequence
        qual = get_ascii_quality(read.query_qualities)

        if read.is_reverse:
            seq = get_complement(seq)[::-1]
            qual = qual[::-1]

        if not read.is_read1 and not read.is_read2:

            if unpaired_read_count >= num_reads:
                chunk += 1
                unpaired_string = 'unpaired_{}.fastq.gz'.format(chunk)
                unpaired_process = spawn_gzip(os.path.join(output_dir, unpaired_string))
                unpaired_read_count = 0
                chunk_files.append(unpaired_string)

            unpaired_process.stdin.write('@{0}\n{1}\n+\n{2}\n'.format(qname, seq, qual))

            unpaired_read_count += 1

            continue

        if qname not in paired_dict:
            paired_dict[qname] = [ None, None ]

        if read.is_read1:
            paired_dict[qname][0] = (seq, qual)

        elif read.is_read2:
            paired_dict[qname][1] = (seq, qual)

        if all(paired_dict[qname]):

            if read_count >= num_reads:
                chunk += 1
                paired_string = 'paired_{}.fastq.gz'.format(chunk)
                paired_process = spawn_gzip(os.path.join(output_dir, paired_string))
                read_count = 0
                chunk_files.append(paired_string)

            v = [qname, paired_dict[qname][0][0], paired_dict[qname][0][1],
                 paired_dict[qname][1][0], paired_dict[qname][1][1]]

            p = '@{0}/1\n{1}\n+\n{2}\n@{0}/2\n{3}\n+\n{4}\n'.format(*v)

            paired_process.stdin.write(p)

            read_count += 2

            del paired_dict[qname]

    sam.close()

    if len(paired_dict.keys()):

        read_count = 0

        for qname in paired_dict.keys():

            if read_count >= num_reads:
                chunk += 1
                unpaired_string = 'unpaired_{}.fastq.gz'.format(chunk)
                unpaired_process = spawn_gzip(os.path.join(output_dir, unpaired_string))
                read_count = 0
                chunk_files.append(unpaired_string)

            #t = next(i for i in paired_dict[qname] if i is not None)

            t = None

            for i in paired_dict[qname]:
                if i:
                    t = i
                    break

            v = [qname, t[0], t[1]]

            u = '@{0}\n{1}\n+\n{2}\n'.format(*v)

            unpaired_process.stdin.write(u)

            read_count += 1

    # Write to interval file
    o = ''.join([ str(c)[:-len('.fastq.gz')] + '\n' for c in chunk_files ])

    interval_file.write(o)

    interval_file.close()

    return

def spawn_gzip(filename):
    p = subprocess.Popen('gzip > ' + filename,
                         bufsize=-1,
                         shell=True,
                         stdin=subprocess.PIPE)
    return p

def get_complement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join([seq_dict[base] for base in seq])

def get_ascii_quality(char_array):
    q = [ chr(c + 33) for c in char_array ]
    q = ''.join(q)
    return q

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help='BAM to convert to FASTQ.')
    parser.add_argument('output_dir', default='./',
                        help='Output directory of the FASTQ files.')
    parser.add_argument('interval_file', type=argparse.FileType('w'),
                        help='Interval file that will be created.')
    parser.add_argument('--num_reads', '-n', type=int, default=75000000,
                        help='Maximum number of reads per FASTQ file.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()
