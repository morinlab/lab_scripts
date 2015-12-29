#!/usr/bin/env python

"""
indelUtils.py
=============
Functions relevant to k-mer and alignment methods
for quantifying indel VAF.
Written by Bruno Grande.
"""

from __future__ import division

from collections import defaultdict
import logging

import numpy as np


# Sequence Functions

def get_seqs(ref_file, chrom, pos, ref, alt, margin):
    """Obtain reference and alternate sequences
    from Ensembl.

    Returns (ref_seq, alt_seq) tuple
    """
    # Calculate start and end positions
    start = max(pos - margin, 0)
    end = pos + margin
    # Extract reference sequence
    ref_seq = ref_file.fetch(reference=chrom, start=start, end=end).decode("utf-8")
    # Strip away any gaps when calculating length
    ref_len = len(ref.strip("-"))
    alt_len = len(alt.strip("-"))
    # Categorize the variant
    if ref_len < alt_len:  # Insertion
        prefix = ref_seq[:margin]
        suffix = ref_seq[margin:]
        alt_seq = prefix + alt + suffix
    elif ref_len > alt_len:  # Deletion
        prefix = ref_seq[:margin-1]
        suffix = ref_seq[margin-1+len(ref):]
        alt_seq = prefix + suffix
    else:  # SNP
        prefix = ref_seq[:margin]
        suffix = ref_seq[margin+1:]
        alt_seq = prefix + alt + suffix
    return ref_seq, alt_seq


def rev_comp(seq):
    """Return reverse complement"""
    cbases = {"A": "T",
              "T": "A",
              "G": "C",
              "C": "G",
              "N": "N"}
    comp = ""
    for base in seq[::-1]:
        comp += cbases[base]
    return comp


# K-mer Functions

class SeqIndexSet(object):

    def __init__(self, seq):
        self.seq = seq
        self.kmer_idxs = {}

    def get_idx(self, k, step, ival):
        """Return k-mer index. Create it if not
        precomputed.
        """
        # Create param key
        key = (k, step, ival)
        # Check if precomputed
        if key in self.kmer_idxs:
            idx = self.kmer_idxs[key]
        else:
            # Create a new index
            idx = defaultdict(set)
            for offset, kmer in kmer_iter(self.seq, k, step, ival):
                idx[kmer].add(offset)
            # Convert to regular dict
            # This will raise KeyErrors if ever a k-mer that doesn't exist is accessed
            # Instead of silently adding the default value
            idx = dict(idx)
            # Store it for later
            self.kmer_idxs[key] = idx
        return idx


def kmer_iter(seq, k, step, ival):
    """Iterate over k-mers using the same
    subsequence pattern.

    Yields (offset, kmer).
    """
    num_kmers = (len(seq) - k * ival)//step + 1
    gap = "-" * (ival - 1)
    for i in range(num_kmers):
        start = i*step
        end = i*step+k*ival
        kmer = gap.join([seq[i] for i in range(start, end, ival)])
        yield start, kmer


def kmer_count(seq, offset, seq_idxs, ref_seq, k, step, ival, min_olap=2):
    """Returns score for k-mers present
    in the given k-mer index.

    Returns the count/score.
    """
    kmer_count = 0
    seq_idx = seq_idxs.get_idx(k=k, step=1, ival=ival)
    # logging.debug("seq_idx: {}".format(seq_idx))
    for start, kmer in kmer_iter(seq, k, step, ival):
        # If offset is set, check for overlap
        if offset and is_overlap(kmer, ref_seq, offset + start, min_olap=min_olap):
            logging.debug("overlapping kmer: {}".format(kmer))
            if kmer in seq_idx:
                kmer_count += 1
        # If offset is not set, simply check if kmer in idx
        elif not offset and kmer in seq_idx:
            kmer_count += 1
    return kmer_count


def calc_kmer_delta(read_seq, offset, ref_idxs, alt_idxs, k, min_delta=1, max_ival=3, min_olap=2):
    """Determines whether read has more k-mers
    in common with reference sequence or alternate
    sequence.

    abs(difference) >= min_delta
    Attempts with interval lengths <= max_ival

    Returns delta in score between the two.
    If positive, aligns better to reference.
    If negative, aligns better to alternate.
    If zero, abs(difference) < min_delta
    """
    ival = 1
    ref_score = 0
    alt_score = 0
    ref_seq = ref_idxs.seq
    logging.debug("calculating kmer delta...")
    while (abs(ref_score - alt_score) < min_delta) and ival <= max_ival:
        logging.debug("ival: {}".format(ival))
        # Find ref scores for forward and reverse and take max
        logging.debug("calculating score for ref...")
        ref_score += kmer_count(read_seq, offset, ref_idxs, ref_seq, k=k, step=1, ival=ival, min_olap=min_olap)
        logging.debug("kmer ref_score: {}".format(ref_score))
        # Find alt scores for forward and reverse and take max
        logging.debug("calculating score for alt...")
        alt_score += kmer_count(read_seq, offset, alt_idxs, ref_seq, k=k, step=1, ival=ival, min_olap=min_olap)
        logging.debug("kmer alt_score: {}".format(alt_score))
        # Increment ival
        ival += 1
    if abs(ref_score - alt_score) < min_delta:
        logging.debug("kmer delta less than min_delta")
        delta = 0
    elif ref_score == 0 and alt_score == 0:
        logging.debug("kmers didn't match anything")
        delta = 0
    else:
        delta = ref_score - alt_score
    logging.debug("delta: {}".format(delta))
    return delta


def is_forward(read_seq, ref_idxs, k, ival):
    """Returns whether read is forward."""
    fread = read_seq
    rread = rev_comp(read_seq)
    fscore = kmer_count(fread, None, ref_idxs, None, k=k, step=1, ival=ival)
    rscore = kmer_count(rread, None, ref_idxs, None, k=k, step=1, ival=ival)
    return fscore > rscore


# Alignment Functions

alphabet = ['A', 'C', 'G', 'T']
score = [[0, 4, 2, 4, 8],
         [4, 0, 4, 2, 8],
         [2, 4, 0, 4, 8],
         [4, 2, 4, 0, 8],
         [8, 8, 8, 8, 8]]


def aln_score(read, ref, offset=None, margin=5):

    # Normalize input strings according to alphabet
    read = str(read).upper()
    ref = str(ref).upper()

    # Edit ref if offset is given
    if offset:
        ref = ref[offset-margin:offset+len(read)+margin]

    # Create distance matrix
    D = np.zeros((len(read)+1, len(ref)+1), dtype=np.int)

    # Initialize first row
    for i in range(1, len(ref)+1):
        D[0, i] = 0

    # Initialize first column
    for i in range(1, len(read)+1):
        D[i, 0] = D[i-1, 0] + score[alphabet.index(read[i-1])][-1]

    # Fill rest of the matrix
    for i in range(1, len(read)+1):
        for j in range(1, len(ref)+1):
            distHor = D[i, j-1] + score[-1][alphabet.index(ref[j-1])]
            distVer = D[i-1, j] + score[alphabet.index(read[i-1])][-1]
            distDiag = D[i-1, j-1] + score[alphabet.index(read[i-1])][alphabet.index(ref[j-1])]
            D[i][j] = min(distHor, distVer, distDiag)

    # Return min of bottom row
    return min(D[-1])


def calc_aln_delta(read_seq, ref_seq, alt_seq, min_delta=8, offset=None, margin=None):
    """Calculate difference in score between
    a local alignment to the reference sequence
    and one to the alternate sequence.

    Returns the difference in score.
    If positive, alignment to reference is better
    If negative, alignment to alternate is better
    If zero, abs(difference) < min_delta
    """
    logging.debug("calculating aln delta...")
    ref_score = aln_score(read_seq, ref_seq, offset=offset, margin=margin)
    logging.debug("aln ref_score: {}".format(ref_score))
    alt_score = aln_score(read_seq, alt_seq, offset=offset, margin=margin)
    logging.debug("aln alt_score: {}".format(alt_score))
    if abs(ref_score - alt_score) < min_delta:
        logging.debug("aln delta less than min_delta")
        delta = 0
    else:
        delta = -(ref_score - alt_score)
    return delta


# Seeding Functions

def find_offset(read, ref_idxs, indel_len, k, max_ival=3, min_delta=5):
    """Find offset of pattern p in k-mer index.

    Returns offset as int.
    """
    ival = 1
    step = 1
    offset_support = defaultdict(int)
    while ival <= max_ival:
        logging.debug("ival: {}".format(ival))
        ref_idx = ref_idxs.get_idx(k, step, ival)
        # Calculate all offset support for this ival
        for pos, kmer in kmer_iter(read, k, step, ival):
            # Check if kmer is in idx first
            if kmer in ref_idx:
                offsets = ref_idx[kmer]
                for offset in offsets:
                    offset_support[offset - pos] += 1
                    offset_support[offset - pos + indel_len] += 1
        # Check if one stands out
        offsets_sorted = sorted([(support, offset) for offset, support in offset_support.items()], reverse=True)
        logging.debug("offsets_sorted: {}".format(offsets_sorted))
        if len(offsets_sorted) > 1:
            # Check if there is enough difference in support between the top 2
            if offsets_sorted[0][0] - offsets_sorted[1][0] >= min_delta:
                return offsets_sorted[0][1]
            # Check if the offsets are within indel_len of each other
            # This implies that the difference between the top 2 isn't big enough
            elif abs(offsets_sorted[0][1] - offsets_sorted[1][1]) == abs(indel_len):
                # If they are equal, choose accordingly
                if offsets_sorted[0][0] == offsets_sorted[1][0]:
                    if indel_len > 0:
                        return min(offsets_sorted[0][1], offsets_sorted[1][1])
                    elif indel_len < 0:
                        return max(offsets_sorted[0][1], offsets_sorted[1][1])
                # If not, choose the best one
                else:
                    return offsets_sorted[0][1]
        elif len(offsets_sorted) == 1 and offsets_sorted[0][0] >= min_delta:
            return offsets_sorted[0][1]
        # If none stands out, increment ival and re-enter loop
        ival += 1
    return None


def is_overlap(read, ref_seq, offset, min_olap=2):
    """Returns whether read overlaps with
    mutation position.
    """
    mid = len(ref_seq) / 2
    return (offset + min_olap <= mid) and (offset + len(read) - min_olap >= mid)


# Read Iteration Function

def get_olap_reads(reads, ref_idxs, indel_len, k, ival, min_olap):
    """Returns a generator that yields reads that overlap
    the mutation with some processing:
        - Replace Ns with As
        - Reverse read if applicable
    """
    for read in reads:
        logging.debug("")
        logging.debug("read: {}".format(read))
        # Replace Ns with As (workaround)
        read = read.rstrip("\n").replace("N", "A")
        # Reverse read if applicable
        logging.debug("determining read orientation...")
        if not is_forward(read, ref_idxs, k=k, ival=ival):
            logging.debug("read was reversed")
            read = rev_comp(read)
        # Find offset
        logging.debug("determining read offset...")
        offset = find_offset(read, ref_idxs, indel_len, k=max(k//2, 5))
        logging.debug("offset: {}".format(offset))
        if not offset:
            logging.debug("read has no offset; skipping")
            continue
        # Determine if overlaps with mutation position
        logging.debug("determining if read overlaps mutation...")
        if not is_overlap(read, ref_idxs.seq, offset, min_olap=min_olap):
            logging.debug("read does not overlap mutation")
            continue
        # Yield overlapping read and its offset
        yield read, offset


# Orphan k-mer function

def get_orphan_kmers(reads, ref_idxs, alt_idxs, k, min_olap, max_ival=3):
    """Returns a dict of orphan k-mers and their counts.
    """
    ival = 1
    orphans = defaultdict(int)
    while ival <= max_ival:
        # Obtain reference and alternate kmers (ref_kmers)
        ref_idx = ref_idxs.get_idx(k=k, step=1, ival=ival)
        alt_idx = alt_idxs.get_idx(k=k, step=1, ival=ival)
        ref_kmers = ref_idx.viewkeys() | alt_idx.viewkeys()
        for read, offset in get_olap_reads(reads, ref_idxs, k, ival, min_olap):
            # Obtain kmers from read that overlap with mutation
            for start, kmer in kmer_iter(read, k=k, step=1, ival=ival):
                if is_overlap(kmer, ref_idxs.seq, offset + start, min_olap=min_olap):
                    if kmer not in ref_kmers:
                        logging.debug("orphan kmer (start: {}): {}".format(start, kmer))
                        orphans[kmer] += 1
        ival += 1
    return orphans
