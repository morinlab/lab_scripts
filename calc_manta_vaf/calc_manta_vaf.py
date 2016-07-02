#!/usr/bin/env python

"""
calc_manta_vaf.py
=================
This script processes Manta VCF files by adding tumour and normal VAF
information to the INFO and sample fields.

Inputs
------
* Manta VCF file

Outputs
-------
* VAF-annotated VCF file

Caveats
-------
* For some reason, the gzip compression done by Manta is incompatible with
  the gzip module (used by pyvcf). Apparently, the magic numbers at the
  beginning of the file don't match what's expected for gzip files. So, for
  the time being, only uncompressed VCF files are supported.
  Tip: you can use "-" as the input VCF file and pipe the uncompressed VCF
  file using zcat through standard input.
"""

from __future__ import print_function, division
import sys
import argparse
from itertools import imap
import vcf as pyvcf


def main():
    args = parse_args()
    vcf = process_vcf(args.vcf)
    vcf_writer = pyvcf.Writer(args.output, template=vcf)
    for record in imap(process_record, vcf):
        vcf_writer.write_record(record)


def parse_args():
    """
    Parse and validate command-line arguments.
    """
    # Command-line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", type=argparse.FileType("r"))
    parser.add_argument("--output", "-o", type=argparse.FileType("w"), default=sys.stdout)
    args = parser.parse_args()
    # No validation required so far
    return args


def process_vcf(vcf_file):
    """
    Convert file object into VCF Reader and edit metadata to allow for
    VAF FORMAT field.
    """
    is_compressed = vcf_file.name.endswith(".gz")
    vcf_reader = pyvcf.Reader(vcf_file, compressed=is_compressed)
    vcf_reader.formats["VAF"] = pyvcf.parser._Format("VAF", 1, "Float", "Variant allele fraction")
    vcf_reader.infos["TVAF"] = pyvcf.parser._Format("TVAF", 1, "Float", "Tumour variant allele fraction")
    vcf_reader.infos["NVAF"] = pyvcf.parser._Format("NVAF", 1, "Float", "Normal variant allele fraction")
    return vcf_reader


def add_vaf(record, tumour_vaf, normal_vaf, ndigits=2):
    """
    Add VAF field to tumour and normal samples in VCF record
    """
    # Round VAF
    tumour_vaf = round(tumour_vaf, ndigits)
    normal_vaf = round(normal_vaf, ndigits)
    # Edit FORMAT and sample fields
    # Create new calldata class that includes VAF
    new_fields = record.FORMAT.split(":") + ["VAF"]
    calldata = pyvcf.model.make_calldata_tuple(new_fields)
    # Update data for tumour and normal
    new_tumour_data = vars(record.samples[1].data)
    new_tumour_data["VAF"] = tumour_vaf
    new_normal_data = vars(record.samples[0].data)
    new_normal_data["VAF"] = normal_vaf
    # Assign new data to samples
    record.samples[1].data = calldata(**new_tumour_data)
    record.samples[0].data = calldata(**new_normal_data)
    # Edit INFO fields
    record.INFO["TVAF"] = tumour_vaf
    record.INFO["NVAF"] = normal_vaf
    return record


def process_record(record):
    """
    Calculate and add VAF information for each sample in the FORMAT
    field.
    """
    TUMOUR_IDX, NORMAL_IDX = 1, 0
    REF_IDX, ALT_IDX = 0, 1
    fields = record.FORMAT.split(":")
    tumour_data = record.samples[TUMOUR_IDX].data
    normal_data = record.samples[NORMAL_IDX].data
    tref, talt, nref, nalt = 0, 0, 0, 0
    if "PR" in fields:
        tref += tumour_data.PR[REF_IDX]
        talt += tumour_data.PR[ALT_IDX]
        nref += normal_data.PR[REF_IDX]
        nalt += normal_data.PR[ALT_IDX]
    if "SR" in fields:
        tref += tumour_data.SR[REF_IDX]
        talt += tumour_data.SR[ALT_IDX]
        nref += normal_data.SR[REF_IDX]
        nalt += normal_data.SR[ALT_IDX]
    tvaf = talt / (talt + tref) if talt + tref != 0 else 0
    nvaf = nalt / (nalt + nref) if nalt + nref != 0 else 0
    record = add_vaf(record, tvaf, nvaf)
    return record


if __name__ == '__main__':
    main()
