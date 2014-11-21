#!/usr/bin/python


import pysam
import argparse
import os.path


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_bam', '-i', dest='in_bam',
                        help='Input BAM file')
    parser.add_argument('--out_bam', '-o', dest='out_bam',
                        help='Output BAM where filtered reads will be written')
    parser.add_argument('--chromosome', '-c', dest='chrom',
                        help="chromosome mutation is on")
    parser.add_argument('--position', '-p', dest='position',
                        help='position in genome')
    parser.add_argument('--ref_allele', '-r', dest='reference_allele')
    parser.add_argument('--mut_allele', '-m', dest='mutant_allele')
    args = parser.parse_args()
    in_bam_name, in_bam_ext = os.path.splitext(args.in_bam)
    if 'out_bam' not in args:
        args.out_bam = in_bam_name + '.db_strand.' + in_bam_ext
    print args

    # Open the BAM file to be used for extracting read data
    bamfile = pysam.Samfile(args.in_bam, "rb")

    # Keep track of how many variants we have removed
    rm_counter = 0
    pos = int(args.position)

    # store the base in the aligned position of each F and R read to allow
    # concordant F/R reads to be counted
    read_bases = {}

    for pcol in bamfile.pileup(reference=args.chrom, start=pos, end=pos + 1, stepper='all',
                               max_depth=10000000):
        # Look at the pileup data for the specific position we are interested in
        genome_coord = pcol.pos + 1
        if genome_coord != pos:
            continue
        var_reads = []
        for pread in pcol.pileups:
            rname = pread.alignment.qname
            base = pread.alignment.seq[pread.qpos]
            read1 = pread.alignment.is_read1
            if not rname in read_bases:
                read_bases[rname] = [0, 0]
            if read1:
                read_bases[rname][0] = base
            else:
                read_bases[rname][1] = base

    base_count = {}
    good_pairs = {}

    # print read_bases
    for read in read_bases:
        if read_bases[read][0] == read_bases[read][1]:
            # print "match: %s" % read
            good_pairs[read] = 1
        else:
            continue
        if not read_bases[read][0] in base_count:
            base_count[read_bases[read][0]] = 1
        else:
            base_count[read_bases[read][0]] += 1

    print base_count

    bamfile_out = pysam.Samfile(args.out_bam, "wb", template=bamfile)

    for read in bamfile.fetch(args.chrom, pos, pos + 1):
        if read.qname in good_pairs:
            bamfile_out.write(read)

    allele_frac = (float(base_count[args.mutant_allele]) /
                   (float(base_count[args.mutant_allele]) +
                    float(base_count[args.reference_allele])))
    print "%s %s %1.9f" % (base_count[args.mutant_allele], base_count[args.reference_allele],
                           allele_frac)


if __name__ == '__main__':
    main()
