import argparse
import os
import pysam
import glob
import gzip as gzip_module
import subprocess

def main():
    args = parse_args()
    og_bam = args.original_bam
    read_set_out = args.read_set_out
    read_set = args.read_set
    files_to_delete = args.files_to_delete
    files_to_delete_path = args.files_to_delete_path
    new_file_path = args.new_file_path
    new_bams = args.new_bam_filename
    new_fastqs = args.new_fastq_filename

    print dir(gzip_module)

    bam_files = []
    fastq_files = []

    original_reads = set()
    new_reads = set()

    to_remove = []

    if og_bam == None and read_set == None:
        raise ValueError('--original_bam or --read_set parameter must be set.')

    if not og_bam == None and read_set_out == None:
        raise ValueError('--read_set_out parameter must be set when using --original_bam.')

    if og_bam == None and not read_set_out == None:
        raise ValueError('--original_bam parameter must be set when using --read_set_out.')

    if not og_bam == None and not read_set == None:
        raise ValueError('--original_bam and --read_set parameters cannot both be set. Only use one.')

    if not new_bams == None and not new_fastqs == None:
        raise ValueError('--new_bam and --new_fastq parameters cannot both be set. Only use one.')

    if new_bams == None and new_fastqs == None:
        raise ValueError('Either --new_bam_filename or --new_fastq_filename must be set.')

    if not files_to_delete == None and files_to_delete_path == None:
        print "Using './' path as default filepath for files indicated by --files_to_delete parameter."

    if new_file_path == None:
        print "Using './' path as default filepath for new BAMs or new FASTQs."

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

        p = subprocess.Popen('gzip > ' + read_set_out + '.gz', bufsize=-1,
                             shell=True, stdin=subprocess.PIPE)
        p.stdin.write('\n'.join(original_reads))

    elif og_bam == None and not read_set == None:
        # read set parameter set and original bam parameter not set.
        if '.gz' in os.path.basename(read_set):
            with gzip_module.open(read_set, 'rb') as f:
                original_reads = set(r.rstrip() for r in f)
        else:
            with open(read_set) as f:
                original_reads = set(r.rstrip() for r in f)

    # Parse new bam or new fastqs into set
    if not new_bams == None and new_fastqs == None:
        if len(new_bams) == 1:
            if '*' in new_bams[0]:
                bam_files = glob.glob(os.path.join(new_file_path, os.path.basename(new_bams[0])))
            else:
                bam_files.append(os.path.join(new_file_path, os.path.basename(new_bams[0])))
        else:
            bam_files = [os.path.join(new_file_path, os.path.basename(b)) for b in new_bams]

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
                fastq_files = glob.glob(os.path.join(new_file_path, os.path.basename(new_fastqs[0])))
            else:
                fastq_files.append(os.path.join(new_file_path, os.path.basename(new_fastqs[0])))
        else:
            fastq_files = [os.path.join(new_file_path, os.path.basename(f)) for f in new_fastqs]

        for fastq_file in fastq_files:
            abs_fastq_file = os.path.abspath(fastq_file)
            if '.gz' in fastq_file:
                with gzip_module.open(abs_fastq_file) as fastq:
                    i = 0
                    for line in fastq:
                        if i % 4 == 0:
                            qname = line.rstrip()
                            new_reads.add(qname[1:])
                        i += 1
            else:
                with open(abs_fastq_file) as fastq:
                    i = 0
                    for line in fastq:
                        if i % 4 == 0:
                            qname = line.rstrip()
                            new_reads.add(qname[1:])
                        i += 1

    # Check if the two sets are the same.
    if not len(original_reads - new_reads) == 0:
        print original_reads - new_reads
        raise ValueError("FAIL: Read names are missing in the new files compared to original file.")
    elif not len(new_reads - original_reads) == 0:
        print new_reads - original_reads
        raise ValueError("FAIL: More read names in new file(s) compared to original file.")
    else:
        print "SUCCESS: Read names in original file matches read names in new file(s)."
        to_remove = []
        if not files_to_delete == None:
            if len(files_to_delete) == 1:
                if '*' in files_to_delete[0]:
                    to_remove = glob.glob(os.path.join(files_to_delete_path, os.path.basename(files_to_delete)))
                else:
                    to_remove.append(os.path.join(files_to_delete_path, os.path.basename(files_to_delete[0])))
            else:
                to_remove = [ os.path.join(files_to_delete_path, os.path.basename(f)) for f in files_to_delete ]
            subprocess.call(["rm"] + to_remove)

    return


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--original_bam', nargs='?',
                        help='BAM that will be realigned.')
    parser.add_argument('--read_set_out', nargs='?',
                        help='Filename of serialized read name set of original bam file.')
    parser.add_argument('--read_set', nargs='?',
                        help='File containing reads of the original BAM.')
    parser.add_argument('--files_to_delete_path', nargs='?', default='./',
                        help='Path to directory of files to delete indicated by --files_to_delete \
                        parameter. [./]')
    parser.add_argument('--files_to_delete', nargs='*',
                        help='Filename(s) or glob of files to delete.')
    parser.add_argument('--new_file_path', nargs='?', default = './',
                        help='Path to directory of location of new files (BAMs or FASTQs). [./]')
    parser.add_argument('--new_bam_filename', nargs='*',
                        help='BAM filename(s) or glob that will undergo read counting to \
                        check for read integrity.')
    parser.add_argument('--new_fastq_filename', nargs='*',
                        help='FASTQ filename(s) or glob that will undergo read counting to \
                        check for read integrity.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
