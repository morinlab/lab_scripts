
import pysam
import argparse
import hashlib

# Read command line arguments
parser = argparse.ArgumentParser()

parser.add_argument('--in_bam', '-i', dest='in_bam', help='Input BAM file')
parser.add_argument('--out_bam','-o', dest='out_bam', help='Output BAM where filtered reads will be written')
parser.add_argument('--chromosome','-c',dest='chrom',help="chromosome mutation is on")
parser.add_argument('--position','-p',dest='position',help='position in genome')
parser.add_argument('--ref_allele','-r',dest='reference_allele')
parser.add_argument('--mut_allele','-m',dest='mutant_allele')

args = parser.parse_args()

print args

# Open the BAM file to be used for extracting read data
bamfile = pysam.Samfile(args.in_bam, "rb")

# Keep track of how many variants we have removed
rm_counter = 0
pos = int(args.position)
read_bases = {} #store the base in the aligned position of each F and R read to allow concordant F/R reads to be counted




for pcol in bamfile.pileup(reference=args.chrom, start=pos, end=pos+1, stepper = 'all', max_depth = 10000000):
	# Look at the pileup data for the specific position we are interested in
	genome_coord = pcol.pos + 1
	if genome_coord != pos:
		continue
	var_reads = []
	for pread in pcol.pileups:
		rname = pread.alignment.qname
		base = pread.alignment.seq[pread.qpos]
		read1 = pread.alignment.is_read1
		if not read_bases.has_key(rname):
			read_bases[rname] = [0,0]
		if read1:
			read_bases[rname][0] = base
		else:
			read_bases[rname][1] = base


base_count = {}

good_pairs = {}
#print read_bases
for read in read_bases:
	if read_bases[read][0] == read_bases[read][1]:
		#print "match: %s" % read
		good_pairs[read] = 1
	else:
		continue
	if not base_count.has_key(read_bases[read][0]):
		base_count[read_bases[read][0]]=1
	else:
		base_count[read_bases[read][0]]+=1

print base_count

bamfile_out = pysam.Samfile(args.out_bam, "wb", template=bamfile)

for read in bamfile.fetch(args.chrom,pos,pos+1):
	if good_pairs.has_key(read.qname):
		bamfile_out.write(read)

allele_frac =  float(base_count[args.mutant_allele])/(float(base_count[args.mutant_allele])+float(base_count[args.reference_allele]))
print "%s %s %1.9f" % (base_count[args.mutant_allele], base_count[args.reference_allele], allele_frac)
