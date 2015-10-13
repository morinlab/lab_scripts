import argparse
import os
<<<<<<< HEAD
import subprocess
import sys

def main():

=======

def main():
    
>>>>>>> 38fd2ab825c276a01774bc3b09bff30ff580b649
    args = parse_args()
    chnk_size = args.num_reads
    bam_infile = args.bam[0]
    paired_outdir = args.paired_outdir[0]
    unpaired_outdir = args.unpaired_outdir[0]

    p_outstring = '@{0}/1\n{1}\n+\n{2}\n@{0}/2\n{3}\n+\n{4}\n'
    up_outstring = '@{0}\n{1}\n+\n{2}\n'

    p_chnk_num = 1
    p_read_count = 0
    up_chnk_num = 1
    up_read_count = 0

    bits = [16, 64, 128, 2048]
<<<<<<< HEAD
=======

>>>>>>> 38fd2ab825c276a01774bc3b09bff30ff580b649
    if not os.path.exists(paired_outdir):
        os.mkdir(paired_outdir)

    if not os.path.exists(unpaired_outdir):
        os.mkdir(unpaired_outdir)

<<<<<<< HEAD
    gzip_cmd = "gzip > "
    p_sp = spawn_gzip(p_chnk_num, paired_outdir, gzip_cmd)
    up_sp = spawn_gzip(up_chnk_num, unpaired_outdir, gzip_cmd)
=======
    p_chnk_file = open_chunk_file(p_chnk_num, paired_outdir)
    up_chnk_file = open_chunk_file(up_chnk_num, unpaired_outdir)

>>>>>>> 38fd2ab825c276a01774bc3b09bff30ff580b649
    # prev_line = [qname, seq, quality, bit 16, bit 64, bit 128]
    prev_line = []
    for line in bam_infile:
        bam_line = line.split('\t')
        sam_flag_list = decompose_flag(bam_line[1], bits)

        # Supplementary alignment check
        if sam_flag_list[3]:
            continue

        if sam_flag_list[0]:
            bam_line[9] = revcomp(bam_line[9])
<<<<<<< HEAD
            bam_line[10] = str(bam_line[10])[::-1]
=======
            bam_line[10] =str( bam_line[10])[::-1]
>>>>>>> 38fd2ab825c276a01774bc3b09bff30ff580b649

        if not len(prev_line):
            prev_line = [bam_line[0]] + bam_line[9:11] + sam_flag_list[0:3]
            continue

        if str(bam_line[0]) == str(prev_line[0]):
<<<<<<< HEAD
            # check readcount, get process, write to its stdin
            t = check_read_count(p_read_count, p_chnk_num, paired_outdir, chnk_size, p_sp, gzip_cmd)
            p_read_count, p_chnk_num, p_sp = t
            p_read_count += 2
            val = []
            if sam_flag_list[1]:
                # curr read is first in pair
                val = [bam_line[0], bam_line[9], bam_line[10],\
                       prev_line[1], prev_line[2]]
            else:
                # curr read is second in pair
                val = [prev_line[0], prev_line[1], prev_line[2],\
                       bam_line[9], bam_line[10]]
            p_sp.stdin.write(p_outstring.format(*val))
            prev_line = []
        else:
            # save prev_line as unpaired read, curr line becomes prev_line
            t = check_read_count(up_read_count, up_chnk_num, unpaired_outdir, chnk_size, up_sp, gzip_cmd)
            up_read_count, up_chnk_num, up_sp = t
            up_read_count += 1
            val = [prev_line[0], prev_line[1], prev_line[2]]
            up_sp.stdin.write(up_outstring.format(*val))
=======
            t = check_read_count(p_read_count, p_chnk_num, p_chnk_file, paired_outdir, chnk_size)
            p_read_count, p_chnk_num, p_chnk_file = t
            if sam_flag_list[1]:
                # curr read is first in pair
                val = [bam_line[0], bam_line[9], bam_line[10], prev_line[1], prev_line[2]]
                p_chnk_file.write(p_outstring.format(*val))
            else:
                # curr read is second in pair
                val = [prev_line[0], prev_line[1], prev_line[2], bam_line[9], bam_line[10]]
                p_chnk_file.write(p_outstring.format(*val))
            prev_line = []
        else:
            # output prev_line as unpaired read, curr line becomes prev_line
            t = check_read_count(up_read_count, up_chnk_num, up_chnk_file, unpaired_outdir, chnk_size)
            up_read_count, up_chnk_num, up_chnk_file = t
            val = [prev_line[0], prev_line[1], prev_line[2]]
            up_chnk_file.write(up_outstring.format(*val))
>>>>>>> 38fd2ab825c276a01774bc3b09bff30ff580b649
            prev_line = [bam_line[0]] + bam_line[9:11] + sam_flag_list[0:3]
    else:
        bam_infile.close()
        if len(prev_line):
<<<<<<< HEAD
            t = check_read_count(up_read_count, up_chnk_num, unpaired_outdir, chnk_size, up_sp, gzip_cmd)
            up_read_count, up_chnk_num, up_sp = t
            up_read_count += 1
            val = [prev_line[0], prev_line[1], prev_line[2]]
            up_sp.stdin.write(up_outstring.format(*val))
=======
            t = check_read_count(up_read_count, up_chnk_num, up_chnk_file, unpaired_outdir, chnk_size)
            up_read_count, up_chnk_num, up_chnk_file = t
            val = [prev_line[0], prev_line[1], prev_line[2]]
            up_chnk_file.write(up_outstring.format(*val))
>>>>>>> 38fd2ab825c276a01774bc3b09bff30ff580b649


def decompose_flag(flag, bits):
    flag_list = [bool(int(flag) & x) for x in bits]
    return flag_list


<<<<<<< HEAD
def check_read_count(rc, chnk_num, chnk_dir, chnk_size, sp, cmd):
    if rc >= chnk_size:
        chnk_num += 1
        sp = spawn_gzip(chnk_num, chnk_dir, cmd)
        rc = 0
    return [rc, chnk_num, sp]
=======
def check_read_count(rc, chnk_num, chnk_file, chnk_dir, chnk_size):
    if rc >= chnk_size:
        chnk_file.close()
        chnk_num += 1
        rc = 0
        chnk_file = open_chunk_file(chnk_num, chnk_dir)
    return [rc, chnk_num, chnk_file]

        
def open_chunk_file(chnk_num, chnk_dir='.'):
    chnk_file = os.path.join(chnk_dir, str(chnk_num) + '.fastq')
    return open(chnk_file, 'w')
>>>>>>> 38fd2ab825c276a01774bc3b09bff30ff580b649


def revcomp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rev = "".join(complement.get(base, base) for base in reversed(seq))
    return rev


<<<<<<< HEAD
def spawn_gzip(chnk_num, chnk_dir, gzip_cmd):
    s = os.path.join(chnk_dir, str(chnk_num) + '.fastq.gz')
    p = subprocess.Popen(gzip_cmd + s, bufsize=-1, shell=True, stdin=subprocess.PIPE)
    return p


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--num_reads', '-n', type=int, default=5000000000,
=======
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--num_reads', '-n', type=int, default=75000000,
>>>>>>> 38fd2ab825c276a01774bc3b09bff30ff580b649
                        help='Specify number of reads per FASTQ file.')
    parser.add_argument('bam', nargs=1, type=argparse.FileType('r'),
                        help='Specify SAM file')
    parser.add_argument('paired_outdir', default='.', nargs=1)
    parser.add_argument('unpaired_outdir', default='.', nargs=1)
    return parser.parse_args()


if __name__ == '__main__':
    import timeit
    print(timeit.timeit('main()', setup='from __main__ import main', number=1))
