"""
tabulate_genes.py
=================
This script tabulates the most recurrently altered
genes in a cohort VCF file.
"""

import argparse
import sys
from collections import defaultdict
import vcf
from annotate_vcf_on_cohort import *


def main():
    """Main program"""

    # Argument parsing
    args = parse_args()

    # Setup
    vcf_reader = vcf.Reader(args.input_vcf)
    vep_cols = parse_vep_cols(vcf_reader)

    # Build dict of genes with affected samples
    genes = defaultdict(set)

    # Iterate over VCF file
    for record in vcf_reader:
        # Filter on NUM_SAMPLES
        if args.max_samples and record.INFO["NUM_SAMPLES"] > args.max_samples:
            continue
        # Parse VEP output and select the first and only one
        vep_effect = parse_vep(vep_cols, record, tag="TOP_CSQ")[0]
        # Extract gene ID and symbol
        gid, gsymbol = vep_effect["Gene"], vep_effect["SYMBOL"]
        # Extract affected samples
        samples = set([call.sample for call in record.samples if call.gt_type != 0])
        # Add samples to genes dict; using gid and gsymbol for readability
        genes[(gid, gsymbol)].update(samples)

    # Order genes by number of affected samples
    genes_list = [(gene[0], gene[1], len(samples)) for gene, samples in genes.items()]
    genes_list.sort(key=lambda x: x[2], reverse=True)

    # Output sorted gene list
    header = "\t".join(["gene_id", "gene_symbol", "num_samples"]) + "\n"
    args.output.write(header)
    for gene in genes_list:
        line = "\t".join(map(str, gene)) + "\n"
        args.output.write(line)

    # Cleanup
    args.output.close()


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("input_vcf", type=argparse.FileType("r"))
    parser.add_argument("--output", "-o", default=sys.stdout, type=argparse.FileType('w'))
    parser.add_argument("--max_samples", "-m", type=int, help="Max. number of samples allowed")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
