#!/usr/bin/env python2

import warnings

NORMAL_IDX = 0
TUMOR_IDX  = 1


def generate_count_index(vcf_reader):
    count_idx = {}
    for record in vcf_reader:
        key = get_idx_key(record)
        t_ref_count, t_alt_count, n_ref_count, n_alt_count = calc_counts(record)
        strelka_qual = get_strelka_qual(record)
        if t_ref_count is None or strelka_qual is None:
            continue
        else:
            count_idx[key] = (t_ref_count, t_alt_count, n_ref_count, n_alt_count, strelka_qual)
    return count_idx


def get_idx_key(record):
    """Generate an index key based on the chromosome, position and reference
    allele.
    """
    if record.is_snp:
        key = (record.CHROM, record.POS, record.REF)
    elif record.is_indel:
        # Adjust position and alleles for MAF format
        chrom, pos, ref, alt = fix_indel_alleles(record.CHROM, record.POS,
                                                 record.REF, record.ALT)
        key = (chrom, pos, ref)
    else:
        warnings.warn("Encountered variant that isn't a SNV or an indel.")
    return key


def fix_indel_alleles(chrom, pos, ref, alt):
    """Adjust indel alleles to follow the shorter MAF format.
    You can assume that one of the alleles has length of one.

    Example:
                |       Deletion        |       Insertion
                |     ref      alt      |     ref      alt
        --------|-----------------------|-----------------------
        VCF     |     CT       T        |     G        GTA
        MAF     |     T        -        |     -        TA
    """
    if len(ref) == 1:  # Insertion
        ref = "-"
        alt = alt[1:]
    elif len(alt) == 1:  # Deletion
        pos += 1
        alt = "-"
        ref = ref[1:]
    else:
        warnings.warn("None of the indel alleles was length == 1.")
    return chrom, pos, ref, alt


def calc_counts(record):
    """Extract the number of reads that support the reference and alternate
    SNV alleles in the tumour and normal samples.
    """
    if record.is_snp:
        return calc_snv_counts(record)
    elif record.is_indel:
        return calc_indel_counts(record)
    else:
        warnings.warn("Encountered variant that isn't a SNV or an indel.")
        return None, None, None, None


def calc_snv_counts(record):
    """Extract the number of reads that support the reference and alternate
    SNV alleles in the tumour and normal samples.
    """
    # Create keys to extract read counts
    ref_key = "{}U".format(record.REF)
    alt_key = "{}U".format(record.ALT[0])
    # Extract read counts (tier 2, which includes tier 1)
    t_ref_count = record.samples[TUMOR_IDX][ref_key][1]
    t_alt_count = record.samples[TUMOR_IDX][alt_key][1]
    n_ref_count = record.samples[NORMAL_IDX][ref_key][1]
    n_alt_count = record.samples[NORMAL_IDX][alt_key][1]
    return t_ref_count, t_alt_count, n_ref_count, n_alt_count


def calc_indel_counts(record):
    """Extract the number of reads that support the reference and alternate
    indel alleles in the tumour and normal samples.
    """
    # Extract read counts (tier 2, which includes tier 1)
    t_ref_count = record.samples[TUMOR_IDX]["TAR"][1]
    t_alt_count = record.samples[TUMOR_IDX]["TIR"][1]
    n_ref_count = record.samples[NORMAL_IDX]["TAR"][1]
    n_alt_count = record.samples[NORMAL_IDX]["TIR"][1]
    return t_ref_count, t_alt_count, n_ref_count, n_alt_count


def get_strelka_qual(record):
    """Extract Strelka quality score (QSS_NT or QSI_NT)."""
    if record.is_snp:
        qual = record.INFO["QSS_NT"]
    elif record.is_indel:
        qual = record.INFO["QSI_NT"]
    else:
        warnings.warn("Encountered variant that isn't a SNV or an indel.")
        qual = None
    return qual
