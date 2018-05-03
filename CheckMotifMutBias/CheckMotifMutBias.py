#!/usr/bin/env python

# DESCRIPTION
#  Checks for evidence of mutations at a given motif, as well as a designated position in the motif
#
#  After obtaining the sequences for each region of interest, all 4 possible orientations of the input
#  motif is generated. The locations of the motif in each region are then identified.
#  The number of mutations overlapping the motif are then calculated, as well as the number of mutations
#  overlapping the position of interest
#  The probability of mutation enrichment/derichment is the calculated using a fisher's exact test:
#
#    Probability of mutation enrichment/derichment for motif:
#        Number of mutations overlapping motif, Total lengths of motifs
#        Number of mutations in gene,           Gene length
#
#    Probability of mutation enrichment/derichment for motif:
#        Number of mutations overlapping position in motif, Total number of positions cross all motifs
#        Number of mutations overlapping position gene-wide, Total number of positions across region
#
# AUTHOR
#  Christopher Rushton
#

import argparse
import os
import pyfaidx
import sys
import re
from scipy.stats import binom_test, fisher_exact

def isValidFile(file, parser):
	"""
	Ensures that the specified file-path actually exists
	"""
	if not os.path.exists(file):
		raise parser.error("Unable to locate \"%s\": No such file or directory" % file)
	else:
		return file


def getArgs():
	"""
	Processes command line arguments
	"""
	parser = argparse.ArgumentParser(description="Checks for enrichment of mutations at a given motif/site in a motif")

	parser.add_argument("-m", "--maf_file", metavar="MAF", required=True, type=lambda x: isValidFile(x, parser),
						help="Input MAF file, listing mutations in the regions of interest")
	parser.add_argument("-r", "-f", "--reference", metavar="FASTA", required=True, type=lambda x: isValidFile(x, parser),
						help="Reference genome, in FASTA format")
	parser.add_argument("-t", "--targets", metavar="BED", required=True, type=lambda x: isValidFile(x, parser),
						help="A BED3 file (or better) listing the regions of interest")
	parser.add_argument("-s", "--signature", metavar="AATYG", required=True, help="Mutation signature sequence")
	parser.add_argument("-i", "--index", metavar="1", required=True, type=int, help="1-based index of the position of interest in the mutation signature")
	parser.add_argument("-o", "--output", metavar="TSV", required=False, default=sys.stdout, help="Output file to write summary report. [default: stdout]")
	parser.add_argument("--method", metavar="STATS", choices=["fisher_exact", "binomial_exact"], default="binomial_exact",
						help="Method to use to test for increased/decreased frequency of mutations [default: %(default)s]")
	parser.add_argument("--annotate_maf", metavar="MAF", required=False, help="Output MAF file listing all mutations, and if they overlap the signature or site of interest")
	args = parser.parse_args()

	# Sanity checks
	if args.index <= 0 or args.index > len(args.signature):
		parser.error("\'-i/--index\' must be a 1-based index corresponding to the position in the \'-s/--signature\' of interest")

	return args


def loadGeneSequences(targetFile, refFile):
	"""
	Obtain the DNA sequences of each region provided

	:args targetFile: A BED3 file (or better) listing the regions of interest
	:args refFile: A file-path to a FASTA file corresponding to the reference genome
	"""

	# Load targets
	targetSeq = {}
	geneNames = {}

	# Load the reference genome
	refGenome = pyfaidx.Fasta(refFile)

	with open(targetFile) as f:
			for line in f:
				line = line.rstrip()
				# Obtain genomic coordinates
				try:
					cols = line.split("\t")
					chrom, start, end = cols[0:3]

					# Generate a unique ID which will be used to identify this sequence
					idName = chrom + ":" + start + "-" + end

					start = int(start) # Account for 1-based indexing
					end = int(end)

				except IndexError:  # This record is malformed or missing. Don't process it
					sys.stderr.write("WARNING: BED entry \'" + line + "\' appears malformed. Ignoring...")

				# Parse the gene name (if present)
				try:
					geneName = cols[3]
				except IndexError:
					# Use the position as a replacement ID
					geneName = idName

				# Obtain the DNA sequence for this gene
				geneSeq = refGenome[chrom][start:end].seq

				# Sanity check: If this region already exists (i.e. it was specified in the BED file multiple times)
				# inform the user
				if idName in targetSeq:
					sys.stderr.write("WARNING: BED entry \'" + line + "\' appears to be a duplicate. Ignoring...")
				else:
					geneNames[idName] = geneName
					targetSeq[idName] = geneSeq

	return targetSeq, geneNames


def loadMutations(mafFile):
	"""
	Stores the position of mutations from an input file in a dictionary

	:param mafFile: A string containing a file-path to the input file
	:return: A dictionary storing chrom: [mutation1, mutation2] (sorted)
	"""

	mutations = {}
	samples = set()
	with open(mafFile) as f:
		for line in f:
			# Skip header lines
			if line.startswith("#") or line.startswith("Hugo_Symbol"):
				continue

			# Parse the chromosome and position of each mutation
			try:
				cols = line.split("\t")
				chrom = cols[4]
				position = int(cols[5]) - 1
				sample = cols[15]
			except (IndexError, TypeError):
				raise TypeError("Unable to parse input MAF file. It may be malformed.")
			if chrom not in mutations:
				mutations[chrom] = []
			mutations[chrom].append(position)

			# Add the sample
			if sample not in samples:
				samples.add(sample)

	for chrom, positions in mutations.items():
		mutations[chrom] = sorted(positions)

	return mutations, len(samples)


def compareMutRate(targetSeq, geneNames, signatures, siteRegex, baseIndexes, mutations, sampleNum, outFile, method):
	"""
	:param targetSeq: A dictionary containing geneID: DNA_Sequence
	:param geneNames: A dictionary storing geneID: GeneName  (ex. chr1:100-1000 : ID3)
	:param signatures: An iterable containing all possible orientations of the mutation signatures
	:param siteRegex: A string representing a regular expression which represents the position of interest in the motif
	:param baseIndexes: An iterable listing indexes for a particular base in signatures
	:param mutations: A dictionary listing chromosome: [position1, position2, position3] for all mutations in the input file
	:param sampleNum: An int listing the number of samples in the input
	:param outFile: A string containing a file-path to an output file, or sys.stdout
	:param method: The method used to calculate the mutation bias (ex. fishers_exact)
	"""

	# Compile the motif signatures
	compiledSig = tuple(re.compile(x) for x in signatures)

	# Store mutation which overlap the motif and site of interest
	totalMotifMut = {}
	totalSiteMut = {}

	with open(outFile, "w") if isinstance(outFile, str) else sys.stdout as o:
		# Write the file header
		header = "\t".join(["gene", "motif_bases", "gene_bases", "gene_mutations", "expected_motif_mut",
							"actual_motif_mut", "motif_bias", "motif_odds_ratio", "motif_site_bases", "gene_site_bases",
							"expected_motif_site_mut", "actual_motif_site_mut", "motif_site_bias", "motif_site_odds_ratio"])
		print(header, file=o)
		for gene, geneSeq in targetSeq.items():

			regionLength = len(geneSeq)
			chrom, pos = gene.split(":")
			geneStart, geneEnd = pos.split("-")
			geneStart = int(geneStart)
			geneEnd = int(geneEnd)
			# First thing we need to is to determine how many times the specified signature occurs
			# in each gene
			sigLoc = {}  # Will store the index of a given match, and the index of the base of interest for that hit

			# Find all matches (in all orientations) for this mutation signature in this gene
			for sig, posIndex in zip(compiledSig, baseIndexes):

				# Find all instances of this motif
				pos = 0

				# This implementation is required to find overlapping motifs
				while True:
					hit = sig.search(geneSeq, pos)
					if hit is None:  # No more instances of this motif exist in the gene sequence
						break
					x = hit.start()
					if x not in sigLoc:
						sigLoc[x] = []
					sigLoc[x].append((hit, hit.start() + posIndex))
					pos = x + 1

			# Calculate the number of bases which overlap the motif
			# Also identify which mutations overlap the signatures
			motifBases = 0
			residueBases = 0  # How many bases overlap this residue?
			motifMut = 0
			nonMotifMut = 0
			residueMut = 0
			nonResidueMut = 0

			if chrom not in totalMotifMut:
				totalMotifMut[chrom] = []
				totalSiteMut[chrom] = []

			# What mutation overlap this gene
			mutOnGene = list(x - geneStart for x in mutations[chrom] if geneStart <= x < geneEnd)
			mutNum = len(mutOnGene)
			currentMut = 0
			positions = sigLoc.keys()
			currentPos = 0
			coveredResidues = set()
			for key in sorted(positions):
				hits = sorted(sigLoc[key], key=lambda x: x[0].end())
				for hit in hits:
					# Add the UNIQUE bases which overlap this hit to the motif overlap count
					if hit[0].start() > currentPos:
						start = hit[0].start()
					else:
						start = currentPos
					currentPos = hit[0].end()
					motifBases += hit[0].end() - start

					# Add the UNIQUE bases overlapping the residue of interest to that count
					if hit[1] not in coveredResidues:
						coveredResidues.add(hit[1])
						residueBases += 1

					# Find any mutations which overlap these bases
					while currentMut < mutNum and mutOnGene[currentMut] < start:
						currentMut += 1
						nonMotifMut += 1

					while currentMut < mutNum and mutOnGene[currentMut] < hit[0].end():
						totalMotifMut[chrom].append(mutOnGene[currentMut] + geneStart)
						motifMut += 1
						currentMut += 1

			while currentMut < mutNum:
				nonMotifMut += 1
				currentMut += 1

			# Count up the number of mutations which overlap the residue of interest
			for mut in mutOnGene:
				if mut in coveredResidues:
					residueMut += 1
					totalSiteMut[chrom].append(mut + geneStart)
				else:
					nonResidueMut += 1

			assert residueMut + nonResidueMut == mutNum
			assert motifMut + nonMotifMut == mutNum

			numCoveredResidues = len(coveredResidues)

			# Calculate the number of positions across the entire gene which correspond to the position of interest
			totalSiteBases = set(x.start() for x in re.finditer(siteRegex, geneSeq))
			numTotalSiteBases = len(totalSiteBases)

			# Calculate the number of mutations which overlap the oosition of interest
			residueGeneMut = list(x for x in mutOnGene if x in totalSiteBases)
			numResidueGeneMut = len(residueGeneMut)

			# Calculate the number of mutations which overlap said positions
			expectedSiteMut = round(float(numCoveredResidues) / float(numTotalSiteBases) * numResidueGeneMut)
			# Calculate the expected number of mutations which are expected to overlap the motif
			expectedMotifMut = round(mutNum * float(motifBases) / float(regionLength))

			# Calculate the p-values for each method
			if method == "fisher_exact":
				motifPVal = fisher_exact([[motifMut, motifBases * sampleNum], [mutNum, regionLength * sampleNum]])[1]
				sitePVal = fisher_exact([[residueMut, numCoveredResidues * sampleNum], [numResidueGeneMut, numTotalSiteBases * sampleNum]])[1]
			elif method == "binomial_exact":
				# Calculate the mutation rate per site across the entire region
				expMutRate = float(mutNum) / (float(regionLength) * sampleNum)
				# Run a binomial test to see if mutations are enriched in the motif
				motifPVal = binom_test(motifMut, motifBases * sampleNum, expMutRate)

				# Calculate the mutation rate per position of interest across the entire region
				expPosMutRate = float(numResidueGeneMut) / (float(numTotalSiteBases) * sampleNum)
				sitePVal = binom_test(residueMut, numCoveredResidues * sampleNum, expPosMutRate)

			# Calculate the odds ratio for the motif and site
			if mutNum == 0 or motifBases == 0: # There are no mutations in the region of interest
				motifOddsRatio = "NA"
			else:
				motifOddsRatio = (motifMut / (motifBases * sampleNum - motifMut)) / (float(mutNum) / (regionLength * sampleNum - float(mutNum)))
			if numResidueGeneMut == 0 or numCoveredResidues ==0 :  # There are no residue mutations in the region
				siteOddsRatio = "NA"
			else:
				siteOddsRatio = (residueMut / (numCoveredResidues * sampleNum - residueMut)) / (float(numResidueGeneMut) / (numTotalSiteBases * sampleNum - float(numResidueGeneMut)))

			# Generate output
			outLine = "\t".join([geneNames[gene], str(motifBases), str(regionLength), str(mutNum), str(expectedMotifMut),
								str(motifMut), str(motifPVal), str(motifOddsRatio), str(numCoveredResidues), str(numTotalSiteBases), str(expectedSiteMut),
								str(residueMut), str(sitePVal), str(siteOddsRatio)])

			print(outLine, file=o)

		return totalMotifMut, totalSiteMut


def complimentSignature(sequence, index):
	"""
	Generate all orientations of the specified sequence

	:param sequence: A string consisting of IUPAC bases corresponding to the signature of interest
	:param index: An integer corresponding to a given position in the sequence (1-based)
	:return: A tuple containing all orientations for the input sequence
	"""

	signatures = []
	indexes = []

	compliment = {
		'A': 'T',
		'T': 'A',
		'C': 'G',
		'G': 'C',
		'U': 'A',
		'R': 'Y',
		'Y': 'R',
		'S': 'S',
		'W': 'W',
		'K': 'M',
		'M': 'K',
		'B': 'V',
		'D': 'H',
		'H': 'D',
		'V': 'B',
		'N': 'N'
	}

	IUPACCodeDict = {
		'A': 'A',  # Adenine
		'C': 'C',  # Cytosine
		'G': 'G',  # Guanine
		'T': 'T',  # Thymine
		'R': '[AG]',  # A or G
		'Y': '[CT]',  # C or T
		'S': '[GC]',  # G or C
		'W': '[AT]',  # A or T
		'K': '[GT]',  # G or T
		'M': '[AC]',  # A or C
		'B': '[CGT]',  # C or G or T
		'D': '[AGT]',  # A or G or T
		'H': '[ACT]',  # A or C or T
		'V': '[ACG]',  # A or C or G
		'N': '[ACGT]',  # any base
	}

	# Generate a reverse compliment of this sequence
	try:
		complimentSeq = "".join(list(compliment[x] for x in sequence))
	except KeyError:  # A non-IUPAC base is present in the sequence
		# Find the bad base
		for char in sequence:
			if char not in complimentSeq:
				raise TypeError("Unrecognized IUPAC base: %s" % char)

		# Something has gone horribly wrong!!!!
		raise TypeError("Unrecognized IUPAC base in: %s" % sequence)

	# Convert ambiguous IUPAC bases into regular expressions

	# First orientation: The input orientation
	forwardRegex = list(IUPACCodeDict[x] for x in sequence)
	signatures.append("".join(forwardRegex))
	indexes.append(index - 1)

	# Second orientation: Reverse
	signatures.append("".join(forwardRegex[::-1]))
	indexes.append(len(forwardRegex) - index)

	# Third orientation: Compliment
	reverseRegex = list(IUPACCodeDict[x] for x in complimentSeq)
	signatures.append("".join(reverseRegex))
	indexes.append(index - 1)

	# Fourth orientation: Reverse compliment:
	signatures.append("".join(reverseRegex[::-1]))
	indexes.append(len(reverseRegex) - index)

	# Finally, generate a regular expression just for the position of interest
	sitePos = IUPACCodeDict[sequence[index - 1]]  # Bases which correspond to the IUPAC symbol for the position of interest
	sitePos = sitePos + IUPACCodeDict[compliment[sequence[index - 1]]]  # Reverse compliment
	sitePos = set(sitePos)  # Remove duplicates
	sitePos.discard("[")
	sitePos.discard("]")
	siteRegex = "[" + "".join(list(sitePos)) + "]"
	return signatures, indexes, siteRegex


def writeMutations(motifMut, siteMut, outFile, inFile, signature):
	"""
	Generates one of more MAF files (as specified) which contains mutations overlapping the motif or site of interest

	:param motifMut: A dictionary listing chrom: [positions] for mutations which overlap the signature of interest
	:param motifMut: A dictionary listing chrom: [positions] for mutations which overlap the position of interest
	:param outFile: A string containing a filepath to the output file for annotated mutations
	:param inFile: A string containing a filepath to the original MAF file
	"""

	# Convert the mutations to sets for fast random lookup of mutations
	tmp = {}
	for chrom, variants in motifMut.items():
		tmp[chrom] = set(variants)
	motifMut = tmp
	tmp = {}
	for chrom, variants in siteMut.items():
		tmp[chrom] = set(variants)
	siteMut = tmp

	with open(inFile) as f, open(outFile, "w") as o:

		for line in f:  # Look at each mutation in the input file

			oLine = line.rstrip()
			# Write the header to both output files
			if line.startswith("#") or line.startswith("Hugo_Symbol"):
				oLine = oLine + "\tMutation_Overlap_" + signature + os.linesep
				o.write(oLine)
				continue

			# Parse the chromosome and position of each mutation
			try:
				cols = line.split("\t")
				chrom = cols[4]
				position = int(cols[5]) - 1

			except (IndexError, TypeError):
				raise TypeError("Unable to parse input MAF file. It may be malformed.")

			variantOverlap = "FALSE"
			if chrom in siteMut and position in siteMut[chrom]:
				variantOverlap = "SITE"
			elif chrom in motifMut and position in motifMut[chrom]:
				variantOverlap = "MOTIF"

			oLine = oLine + "\t" + variantOverlap + os.linesep
			o.write(oLine)


def main(args=None):

	if args is None:
		args = getArgs()

	# Load gene sequences
	geneSequences, geneNames = loadGeneSequences(args.targets, args.reference)

	# Generate all orientations of the signature sequence
	signatures, mutIndex, siteRegex = complimentSignature(args.signature, args.index)

	# Load the mutations from the MAF file
	mutations, sampleCount = loadMutations(args.maf_file)

	# Calculates the expected and actual number of mutations that overlap said motif, and the position of interest
	motifMut, siteMut = compareMutRate(geneSequences, geneNames, signatures, siteRegex, mutIndex, mutations, sampleCount, args.output, args.method)

	if args.annotate_maf is not None:
		writeMutations(motifMut, siteMut, args.annotate_maf, args.maf_file, args.signature)


if __name__ == "__main__":
	main()
