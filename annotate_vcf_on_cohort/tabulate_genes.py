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
from annotate_vcf_on_cohort import parse_vep_cols, parse_vep, create_pos_id


def main():
    """Main program"""

    # Argument parsing
    args = parse_args()

    # Setup
    vcf_reader = vcf.Reader(args.input_vcf)
    vep_cols = parse_vep_cols(vcf_reader)

    # Create set of genes to be excluded
    excl_genes_set = build_exclude_genes(args.exclude_genes)

    # Create set of positions to be excluded
    excl_pos_set = build_exclude_positions(args.exclude_positions)

    # Build dict of genes with affected samples
    genes = defaultdict(set)

    # Iterate over VCF file
    for record in vcf_reader:
        # Filter on position, if applicable
        pos_id = create_pos_id(record.CHROM, record.POS)
        if pos_id in excl_pos_set:
            continue
        # Filter on NUM_SAMPLES
        if args.max_samples and record.INFO["NUM_SAMPLES"] > args.max_samples:
            continue
        # Parse VEP output and select the first and only one
        vep_effect = parse_vep(vep_cols, record, tag="TOP_CSQ")[0]
        # Exclude on gene ID or symbol
        if vep_effect["Gene"] in excl_genes_set or vep_effect["SYMBOL"] in excl_genes_set:
            continue
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
    parser.add_argument("--exclude_genes", type=argparse.FileType("r"), help="List of gene IDs or symbols to exclude")
    parser.add_argument("--exclude_positions", type=argparse.FileType("r"), help="List of genomic positions to exclude (format: CHROM\\tPOS)")
    args = parser.parse_args()
    return args


def build_exclude_genes(exclude_file):
    """Build set of genes to be excluded"""
    if exclude_file is None:
        return set()
    lines = exclude_file.readlines()
    return set(lines)


def build_exclude_positions(exclude_file):
    """Build dict of positions to be excluded"""
    positions = set()
    if exclude_file is None:
        return positions
    for line in exclude_file:
        cols = line.split("\t")
        chrom, pos = cols[0:2]
        pos_id = create_pos_id(chrom, pos)
        positions.add(pos_id)
    return positions


if __name__ == '__main__':
    main()
