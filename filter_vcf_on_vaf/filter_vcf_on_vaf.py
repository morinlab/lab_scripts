#!/usr/bin/env python3

import argparse
import numpy as np
from cyvcf2 import VCF

# Tumour sample column index per VCF caller
# MuTect doesn't consistently order the samples
TUMOUR_COL = {
    "strelka": 1,
    "mutect": None
}


def main():
    """Main program"""
    args = parse_args()
    variants = list(args.vcf())
    calc_vaf = get_calc_vaf(args.caller, variants)
    args.output.write(args.vcf.raw_header)
    for variant in variants:
        # Skip anything that isn't a SNP/SNV
        if not variant.is_snp:
            continue
        # Calculate VAF and output if VAF is above threshold
        vaf = calc_vaf(variant)
        if vaf >= args.min_vaf:
            args.output.write(str(variant))


def parse_args():
    """Parse and validate command-line arguments"""
    # Parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", type=VCF, metavar="VCF_FILE")
    parser.add_argument("--output", "-o", type=argparse.FileType("w"), default="-")
    parser.add_argument("--min_vaf", "-m", type=float, default=0.15)
    parser.add_argument("--caller", "-c", choices=["strelka", "mutect"])
    args = parser.parse_args()
    # Validating
    if not args.caller:
        args.caller = detect_caller(args.vcf)
    return args


def detect_caller(vcf):
    """Detect VCF caller for VAF parsing"""
    # Check for Strelka
    if "consolidateResults.pl" in vcf.raw_header:
        caller = "strelka"
    elif "ID=MuTect" in vcf.raw_header:
        caller = "mutect"
    else:
        raise ValueError("Cannot detect VCF caller")
    return caller


def get_calc_vaf(caller, variants):
    """Calculate the tumour VAF"""
    if caller == "strelka":
        calc_vaf = calc_tumour_vaf_strelka()
    elif caller == "mutect":
        calc_vaf = calc_tumour_vaf_mutect(variants)
    else:
        raise ValueError(f"Invalid VCF caller: {caller}")
    return calc_vaf


def calc_tumour_vaf_strelka():
    """Calculate tumour VAF from Strelka VCF files"""
    def calc_tumour_vaf_strelka_helper(variant):
        # Extract allele counts
        ref_key = variant.REF + "U"
        alt_key = variant.ALT[0] + "U"
        # Use tier 2 counts (to be as similar as possible to other callers)
        ref_count = variant.format(ref_key)[TUMOUR_COL["strelka"]][1]
        alt_count = variant.format(alt_key)[TUMOUR_COL["strelka"]][1]
        # Calculate VAF
        vaf = alt_count / (alt_count + ref_count)
        return vaf
    return calc_tumour_vaf_strelka_helper


def calc_tumour_vaf_mutect(variants):
    """Calculate tumour VAF from MuTect VCF files"""
    tumour_col = find_tumour_column(variants)
    def calc_tumour_vaf_mutect_helper(variant):
        # Just use MuTect's FA value in the sample columns
        vaf = variant.format("FA")[tumour_col]
        return vaf
    return calc_tumour_vaf_mutect_helper


def find_tumour_column(variants):
    """Identify tumour column based on average VAF (MuTect only)"""
    vafs = [var.format("FA") for var in variants]
    means = np.concatenate(vafs, axis=1).mean(axis=1)
    if means[0] > means[1]:
        tumour_idx = 0
    else:
        tumour_idx = 1
    return tumour_idx


if __name__ == '__main__':
    main()
