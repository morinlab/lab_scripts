#!/usr/bin/env python

"""
calc_delly_vaf.py
=================
This script prepares a matrix of variant allele fractions (VAFs)
for each structural variant (SV) called by DELLY per related tumour
sample. SVs are to be called jointly by DELLY, such that the same
SVs are genotyped together. Otherwise, SVs would need to be matched
up across multiple VCF files and I prefer to let DELLY deal with
this.

Inputs
------
- DELLY VCF file (any SV type)

Outputs
-------
- Tab-delimited file containing the VAF matrix

Known Issues
------------
- None
"""

import argparse
import re
import vcf


# Here are the standard columns for the output file.
# To this, you must add one column per tumour sample,
# which will contain the VAFs
BASE_COLUMNS = ["sv_id", "chr1", "pos1", "chr2", "pos2", "strand", "type"]


def main():

    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("input_vcf", type=argparse.FileType("r"), help="Input DELLY VCF file")
    parser.add_argument("output", type=argparse.FileType("w"), help="Output matrix TSV file")
    parser.add_argument("--tumour1_regex", "-r1", required=True)
    parser.add_argument("--tumour2_regex", "-r2", required=True)
    parser.add_argument("--min_cov", "-c", type=int, default=0,
                        help="Minimum coverage at SV locus")
    args = parser.parse_args()

    # Argument processing
    t_regex = re.compile(r"({}|{})".format(args.tumour1_regex, args.tumour2_regex))
    invcf = vcf.Reader(args.input_vcf)
    output = args.output
    min_cov = args.min_cov

    # VCF parsing
    is_first = True
    for record in invcf:
        # Create header on first record
        if is_first:
            header = create_header(record, t_regex)
            output.write(header)
            is_first = False
        # Get SV info for first columns
        info_cols = [
            record.ID,
            record.CHROM,
            record.POS,
            record.INFO["CHR2"],
            record.INFO["END"],
            record.INFO["CT"],
            record.INFO["SVTYPE"]]
        # Ensure minimum coverage in all samples
        if all([calc_cov(call) >= min_cov for call in get_tumour_calls(record, t_regex)]):
            # Calculate VAF for each tumour sample
            vaf_cols = []
            for call in get_tumour_calls(record, t_regex):
                vaf = calc_vaf(call)
                vaf_cols.append(vaf)
            # Output
            all_cols = [str(col) for col in info_cols + vaf_cols]
            output.write("\t".join(all_cols) + "\n")

    # Closing files
    args.input_vcf.close()
    output.close()


def create_header(vcf_record, t_regex):
    """Create header for output file based on first VCF record
    """
    tumour_sample_names = [call.sample for call in get_tumour_calls(vcf_record, t_regex)]
    columns = BASE_COLUMNS + tumour_sample_names
    header = "\t".join(columns) + "\n"
    return header


def get_tumour_calls(vcf_record, t_regex):
    """Return calls for tumour samples linked to a given VCF record
    """
    tumour_calls = [call for call in vcf_record.samples if t_regex.search(call.sample)]
    return tumour_calls


def calc_cov(vcf_call):
    """Calculate total coverage for a given DELLY VCF record
    and sample
    """
    cov = (int(vcf_call.data.DR) + int(vcf_call.data.RR) +
           int(vcf_call.data.DV) + int(vcf_call.data.RV))
    return cov


def calc_vaf(vcf_call):
    """Calculate VAF for a given DELLY VCF record
    and sample
    """
    ref_count = float(vcf_call.data.DR) + float(vcf_call.data.RR)
    alt_count = float(vcf_call.data.DV) + float(vcf_call.data.RV)
    if ref_count + alt_count == 0:
        return 0
    vaf = alt_count / (ref_count + alt_count)
    return round(vaf, 4)


if __name__ == "__main__":
    main()
