# useful functions for dealing with BAM files.

import warnings
import subprocess
import tempfile
import sys
import logging
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
#import swalign
import indelUtils


def num_mismatches(ref, query):
    """Number of matches between query and reference."""
    mismatches = 0
    for i in range(len(ref)):
        if query[i] != ref[i]:
            mismatches += 1
    return mismatches


def multi_align(ref1, ref2, queries):
    """
    Do MSA of several query sequences to two references. Return the reference
    index with the fewest mismatches for each sequence.
    """
    tfile = tempfile.NamedTemporaryFile()
    tfile.write(">ref1\n{}\n".format(ref1))
    tfile.write(">ref2\n{}\n".format(ref2))
    n_queries = 0
    for query in queries:
        tfile.write(">query{}\n{}\n".format(n_queries, query))
        n_queries += 1
    tfile.flush()

    outfile = tempfile.TemporaryFile()
    cmd = ["mafft", "--auto", tfile.name]
    subprocess.Popen(cmd, stdout=outfile, stderr=subprocess.PIPE).communicate()
    outfile.seek(0)
    align = SeqIO.parse(outfile, "fasta")

    ref1 = next(align)
    ref2 = next(align)
    scores = [0, 0]
    ties = 0
    for query in align:
        mismatch1 = num_mismatches(ref1, query)
        mismatch2 = num_mismatches(ref2, query)
        if mismatch1 < mismatch2:
            scores[0] += 1
        elif mismatch1 > mismatch2:
            scores[1] += 1
        else:
            ties += 1
    if ties > 0:
        msg = "Discarding {} reads aligning equally to ref and alt".format(ties)
        warnings.warn(msg)

    return scores


def multi_muscle_align(ref1, ref2, queries):
    """
    Do MSA of several query sequences to two references using multiple sequence aligner MUSCLE. Return the reference
    index with the fewest mismatches for each sequence.
    """
    tfile = tempfile.NamedTemporaryFile()
    tfile.write(">ref1\n{}\n".format(ref1))
    tfile.write(">ref2\n{}\n".format(ref2))
    n_queries = 0
    for query in queries:
        tfile.write(">query{}\n{}\n".format(n_queries, query))
        n_queries += 1
    tfile.flush()

    outfile = tempfile.TemporaryFile()
    muscle_cline = MuscleCommandline(tfile.name)
    subprocess.Popen(str(muscle_cline), stdout=outfile, stderr=subprocess.PIPE, shell=(sys.platform != "win32"))
    outfile.seek(0)
    align = SeqIO.parse(outfile, "fasta")
    ref1 = next(align)
    ref2 = next(align)
    scores = [0, 0]
    ties = 0
    for query in align:
        mismatch1 = num_mismatches(ref1, query)
        mismatch2 = num_mismatches(ref2, query)
        if mismatch1 < mismatch2:
            scores[0] += 1
        elif mismatch1 > mismatch2:
            scores[1] += 1
        else:
            ties += 1
    if ties > 0:
        msg = "Discarding {} reads aligning equally to ref and alt".format(ties)
        warnings.warn(msg)

    return scores


def pair_align(ref_seq, alt_seq, reads):
    """
    Does pairwise alignment of read sequences to two references (original and altered). extracts scores for each alignment and counts indel/ref/tie.
    """
    match = 2
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring, full_query=True, gap_penalty=-10, gap_extension_penalty=-0.5)
    n_queries = 0
    scores = [0, 0]
    ties = 0
    for read in reads:

        n_queries += 1
        ref_alignment = sw.align(ref_seq, read)
        alt_alignment = sw.align(alt_seq, read)
        ref_score = ref_alignment.score
        alt_score = alt_alignment.score
        if alt_score > ref_score:
            scores[1] += 1
        elif alt_score < ref_score:
            scores[0] += 1
        else:
            ties += 1
    if ties > 0:
        msg = "Discarding {} reads aligning equally to ref and alt".format(ties)
        warnings.warn(msg)

    return scores


def get_indel_vaf_pileup(samfile, reffile, chromosome, pos, ref, alt, minimum_mapping_qual):
    """This function uses the I/D CIGAR operations at a
    sepcific position to calculate the VAF of an indel
    specified by the user. Only cigar operations that match
    the length of the event will be counted and all other reads
    will be considered equivalent to the 'reference' allele
    I am at this time not able to determine how to match the length
    of an insertion to a pre-specified event so this function currently
    assumes all I and D operations overlapping with a site are equivalent"""
    chromosome_long = "chr" + chromosome
    pad_length = 250

    # first generate a pileup object in the region specified with padding
    region_start = pos-pad_length
    region_end = pos+pad_length
    pos_depths = {}  # store coverage depth at every position in the region
    pos_delcounts = {}  # store number of deletion operations at every position in the region
    pos_inscounts = {}  # store number of insertion operations at every position in the region
    try:
        pileup = samfile.pileup(chromosome_long, region_start, region_end, fastafile=reffile, max_depth=1000000)
    except ValueError:
        pileup = samfile.pileup(chromosome, region_start, region_end, fastafile=reffile, max_depth=1000000)
    for x in pileup:
        genome_coordinate = x.pos+1
        # skip positions in padded region
        if genome_coordinate < region_start or genome_coordinate > region_end:
            continue
        dep = x.n
        if genome_coordinate not in pos_depths:
            pos_depths[genome_coordinate] = dep
        reads = x.pileups
        nref_reads = {"A": {}, "C": {}, "T": {}, "G": {}, "I": {}, "D": {}, "N": {}}
        ref_count = 0
        for read in reads:
            mapq = read.alignment.mapq
            if mapq >= minimum_mapping_qual:
                try:
                    read.alignment.seq[read.query_position]
                except TypeError:
                    pass
                if not read.indel and not read.is_del:
                    continue

                if read.indel > 1:

                    if genome_coordinate not in pos_inscounts:
                        pos_inscounts[genome_coordinate] = 1

                    else:
                        pos_inscounts[genome_coordinate] += 1
                elif read.indel < 1:

                    if genome_coordinate not in pos_delcounts:
                        pos_delcounts[genome_coordinate] = 1

                    else:
                        pos_delcounts[genome_coordinate] += 1

                if read.indel > 1:
                    if read.alignment.qname in nref_reads["I"]:
                        nref_reads["I"][read.alignment.qname] += 1
                    else:
                        nref_reads["I"][read.alignment.qname] = 1
                elif read.is_del:
                    if read.alignment.qname in nref_reads["D"]:
                        nref_reads["D"][read.alignment.qname] += 1
                    else:
                        nref_reads["D"][read.alignment.qname] = 1
                else:
                    ref_count += 1
    nref_count = 0
    for position in range(pos-2, pos+2):
        if position not in pos_depths:
            continue
        if position not in pos_delcounts:
            pos_delcounts[position] = 0
        if position not in pos_inscounts:
            pos_inscounts[position] = 0

        if alt == "-":
            # deletion
            if position == pos:
                nref_count = pos_delcounts[position]
        elif len(alt) > len(ref):
            # insertion
            if position == pos:
                nref_count = pos_inscounts[position]
        elif ref == "-":
            # insertion
            if position == pos:
                nref_count = pos_inscounts[position]
        else:
            logging.error("Is this an indel? %s %s %s" % (pos, ref, alt))
            raise ValueError
    ref_count = pos_depths[pos] - nref_count

    return(ref_count, nref_count)


def kmer_count_and_aln(ref_seq, alt_seq, reads, params={}):
    """This function infers whether a read supports the
    reference or alternate allele and then counts them up.
    It does so by using a hybrid approach combining k-mer
    counting and alignment.

    A first pass is done using k-mer counting. If the read
    has more k-mers in common with the reference than the
    alternate, then it is classified as supporting the
    reference allele, and vice versa.

    If k-mer counting cannot discriminate between the
    reference and alternate alleles, a variant of local
    alignment using dynamic programming is used. First,
    the read is seeded relative to the reference using
    a k-mer index. Second, the read is aligned to that
    portion of the reference sequence. The same is done
    for the alternate sequence. The read is classified
    according to the best alignment score.

    Initial implementation written by Bruno Grande.

    Known Issues
    ------------
    - Bases N aren't handled well. They are simply
      replaced to As.
    """
    # Set parameters
    defaults = {
        "k": 13,
        "ival": 2,
        "min_olap": 0,
        "min_delta_kmer": 3,
        "max_ival": 5,
        "min_delta_aln": 8
    }
    defaults.update(params)
    params = defaults
    logging.debug("params: {}".format(params))
    # Initialize some variables
    ref_count = 0
    alt_count = 0
    amb_count = 0
    indel_len = len(alt_seq) - len(ref_seq)  # + for ins and - for del
    ref_idxs = indelUtils.SeqIndexSet(ref_seq)
    alt_idxs = indelUtils.SeqIndexSet(alt_seq)
    # Log ref and alt sequences
    logging.debug("ref_seq: {}".format(ref_seq))
    logging.debug("alt_seq: {}".format(alt_seq))
    # Iterate over reads
    for read, offset in indelUtils.get_olap_reads(reads, ref_idxs, indel_len, k=params["k"], ival=params["ival"], min_olap=params["min_olap"]):
        # Calculate score delta using k-mer method
        kmer_delta = indelUtils.calc_kmer_delta(read, offset, ref_idxs, alt_idxs, k=params["k"], min_delta=params["min_delta_kmer"], max_ival=params["max_ival"], min_olap=params["min_olap"])
        logging.debug("kmer_delta: {}".format(kmer_delta))
        # if skip_kmer_counting:
        #     kmer_delta = 0
        if kmer_delta > 0:
            ref_count += 1
            logging.debug("read classified as reference by kmer method")
        elif kmer_delta < 0:
            alt_count += 1
            logging.debug("read classified as alternate by kmer method")
        else:
            # If k-mer method can't discriminate between ref and alt, use alignment method
            # Estimate appropriate margin (esp. if insertion)
            margin = max(len(alt_seq) - len(ref_seq) + 5, 5)
            # Calculate score delta using alignment method
            aln_delta = indelUtils.calc_aln_delta(read, ref_seq, alt_seq, min_delta=params["min_delta_aln"], offset=offset, margin=margin)
            logging.debug("aln_delta: {}".format(aln_delta))
            if aln_delta > 0:
                ref_count += 1
                logging.debug("read classified as reference by aln method")
            elif aln_delta < 0:
                alt_count += 1
                logging.debug("read classified as alternate by aln method")
            else:
                amb_count += 1
                logging.debug("read classified as ambiguous")
    return (ref_count, alt_count)


def count_indels(samfile, reffile, chrom, pos, ref, alt, mode, min_mapq=20):
    """
    Count occurences of the reference and alternate indel allele at a given
    position, by alignment score.
    """
    # Log
    logging.debug("")
    logging.debug("")
    logging.debug("")
    logging.debug("indel: {} {} {} {}".format(chrom, pos, ref, alt))

    # Extract reads
    reads = samfile.fetch(chrom, pos, pos+len(ref))
    reads = [r.seq for r in reads if r.mapq >= min_mapq and not r.is_duplicate]

    # If there are no reads, return zero counts
    if len(reads) == 0:
        return {ref: 0, alt: 0}

    # Extract ref and alt sequences
    margin = max(map(len, reads)) + 10  # Dynamically set margin based on read length
    ref_seq, alt_seq = indelUtils.get_seqs(reffile, chrom, pos, ref, alt, margin)

    # Calculate read counts for ref and alt
    method = MODES[mode]
    if mode == "pileup":
        # needs different data using the pileup method
        counts = method(samfile, reffile, chrom, pos, ref, alt, min_mapq)
    else:
        counts = method(ref_seq, alt_seq, reads)
    return {ref: counts[0], alt: counts[1]}


def count_bases_pileup(pileup, position, min_baseq=15, min_mapq=20):
    """
    Count the number of times each nucleotide occurs at a given position in a
    pileup object, subject to minimum base quality and map quality constraints.
    Return a dictionary keyed by nucleotides.
    """
    counts = dict.fromkeys(["A", "T", "C", "G"], 0)
    for x in pileup:
        if position == x.pos + 1:
            for read in x.pileups:
                dup = read.alignment.is_duplicate
                qc_fail = read.alignment.is_qcfail
                low_mapq = read.alignment.mapq < min_mapq
                if not (dup or qc_fail or low_mapq or read.is_del or read.is_refskip):
                    base_qual = ord(read.alignment.qual[read.query_position])-33
                    if base_qual >= min_baseq:
                        base = read.alignment.seq[read.query_position]
                        try:
                            counts[base] += 1
                        except KeyError:
                            pass

    return counts


def count_bases(samfile, reffile, chrom, pos):
    """
    Count the number of times each nucleotide occurs at a given position in a
    bam file. Return a dictionary keyed by nucleotides.
    """
    start = max(pos-200, 0)
    end = pos+200
    reffile.fetch(reference=chrom, start=pos-1, end=pos)
    pileup = samfile.pileup(chrom, start, end, fastafile=reffile)
    return count_bases_pileup(pileup, pos)


MODES = {
    "mafft": multi_align,
    "muscle": multi_muscle_align,
    "swalign": pair_align,
    "hybrid": kmer_count_and_aln,
    "pileup": get_indel_vaf_pileup
}
