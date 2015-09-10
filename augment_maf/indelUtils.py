#!/usr/bin/env python

"""
indelUtils.py
=============
Functions relevant to k-mer and alignment methods
for quantifying indel VAF.
Written by Bruno Grande.
"""

from __future__ import print_function
from __future__ import division

from collections import defaultdict
from random import shuffle

import numpy as np


# Sequence Functions

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
            # Store it for later
            self.kmer_idxs[key] = idx
        return idx


def kmer_iter(seq, k, step, ival):
    """Iterate over k-mers using the same
    subsequence pattern.

    Yields (offset, kmer).
    """
    num_kmers = (len(seq) - k * ival)//step + 1
    kmer_ids = range(num_kmers)
    shuffle(kmer_ids)
    for i in kmer_ids:
        start = i*step
        end = i*step+k*ival
        kmer = seq[start:end:ival]
        yield start, kmer


def kmer_count(seq, kmer_idx, k, step, ival):
    """Returns score for k-mers present
    in the given k-mer index.

    Returns the count/score.
    """
    kmer_count = 0
    for offset, kmer in kmer_iter(seq, k, step, ival):
        if kmer in kmer_idx:
            kmer_count += 1
    return kmer_count


def calc_kmer_delta(read_seq, ref_idxs, alt_idxs, min_delta=1, max_ival=3):
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
    while (abs(ref_score - alt_score) < min_delta) and ival <= max_ival:
        # Generate k-mer indexes for this ival
        ref_idx = ref_idxs.get_idx(k=K, step=1, ival=ival)
        alt_idx = alt_idxs.get_idx(k=K, step=1, ival=ival)
        # Find ref scores for forward and reverse and take max
        ref_score += kmer_count(read_seq, ref_idx, k=K, step=1, ival=ival)
        # Find alt scores for forward and reverse and take max
        alt_score += kmer_count(read_seq, alt_idx, k=K, step=1, ival=ival)
        # Increment ival
        ival += 1
    if abs(ref_score - alt_score) < min_delta:
        delta = 0
    else:
        delta = ref_score - alt_score
    return delta


def is_forward(read_seq, ref_idxs):
    """Returns whether read is forward."""
    fread = read_seq
    rread = rev_comp(read_seq)
    ref_idx = ref_idxs.get_idx(k=K, step=1, ival=2)
    fscore = kmer_count(fread, ref_idx, k=K, step=1, ival=2)
    rscore = kmer_count(rread, ref_idx, k=K, step=1, ival=2)
    return fscore > rscore


# Alignment Functions

alphabet = ['A', 'C', 'G', 'T']
score = [[0, 4, 2, 4, 8],
         [4, 0, 4, 2, 8],
         [2, 4, 0, 4, 8],
         [4, 2, 4, 0, 8],
         [8, 8, 8, 8, 8]]


def aln_score(read, ref, offset=None, margin=5):

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
    ref_score = aln_score(read_seq, ref_seq, offset=offset, margin=margin)
    alt_score = aln_score(read_seq, alt_seq, offset=offset, margin=margin)
    if abs(ref_score - alt_score) < min_delta:
        delta = 0
    else:
        delta = -(ref_score - alt_score)
    return delta


# Seeding Functions

def find_offset(read, ref_idxs, k, step, ival, min_support=3):
    """Find offset of pattern p in k-mer index.

    Returns offset as int.
    """
    offset_support = defaultdict(int)
    ref_idx = ref_idxs.get_idx(k, step, ival)
    for pos, kmer in kmer_iter(read, k, step, ival):
        offsets = ref_idx[kmer]
        for offset in offsets:
            offset_support[offset - pos] += 1
        vals = offset_support.values()
        if any(map(lambda x: x >= min_support, vals)):
            max_support = max(vals)
            best_offsets = [offset for offset, support in offset_support.items() if support == max_support]
            if len(best_offsets) > 1:
                continue
            else:
                return best_offsets[0]
    return None


def is_overlap(read, ref_seq, offset, min_olap=2):
    """Returns whether read overlaps with
    mutation position.
    """
    mid = len(ref_seq) / 2
    return (offset + min_olap <= mid) and (offset + len(read) - min_olap >= mid)
