#!/usr/bin/env python

"""
autoigv_per_gene.py
===================
Create necessary files (shell script and input files)
for autoigv.py to screenshot every position in various
samples, organized per gene.

Requirements
------------
- Python 3

Inputs
------
- VEP-annotated VCF file
- Map sample names to BAM files
"""

import argparse
import os
from collections import defaultdict

import vcf
from annotate_vcf_on_cohort import parse_vep_cols, parse_vep

POSITIONS_NAME = "positions.txt"
SCRIPT_NAME = "run_autoigv.sh"
PREFS_NAME = "autoIGVprefs.ini"


def main():
    args = parse_args()
    bam_map = parse_bam_map_file(args.bam_map_file)
    vcf_reader = vcf.Reader(args.vcf_file)
    vep_cols = parse_vep_cols(vcf_reader)
    gene_list = parse_genes(args.genes)
    records_per_gene = parse_vcf_file(vcf_reader, vep_cols, gene_list)
    for gene, records in records_per_gene.items():
        gene_dir = os.path.join(args.output_dir, gene)
        if not os.path.exists(gene_dir):
            os.mkdir(gene_dir)
        # Create script file
        with open(os.path.join(gene_dir, SCRIPT_NAME), "w") as sf:
            sf.write(generate_autoigv_cmd(args.python, args.autoigv, args.genome))
        # Create positions file
        with open(os.path.join(gene_dir, POSITIONS_NAME), "w") as pf:
            pf.write(generate_autoigv_positions(gene, records, bam_map))
    # Create master script file
    with open(os.path.join(args.output_dir, "run_all.sh"), "w") as mf:
        mf.write(generate_master_script())
    # Create prefs file
    with open(os.path.join(args.output_dir, PREFS_NAME), "w") as pf:
        pf.write(generate_prefs())


def parse_args():
    """ Parse and validate command-line arguments """
    # Parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf_file", type=argparse.FileType("r"), help="Input VCF file")
    parser.add_argument("bam_map_file", type=argparse.FileType("r"), help="Maps sample names (1st col) to BAM files (2nd col)")
    parser.add_argument("--output_dir", "-o", default=".", help="Output directory where all the gene subdirectories will be created")
    parser.add_argument("--autoigv", "-a", default="autoIGV1.5.py", help="Path to autoIGV")
    parser.add_argument("--python", "-p", default="python", help="Path to python")
    parser.add_argument("--genome", "-g", default="hg19", help="Genome to use")
    parser.add_argument("--genes", type=argparse.FileType("r"), help="File containing list of genes to include (Ensembl IDs or symbols)")
    args = parser.parse_args()
    # Validating output directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    # Return
    return args


def parse_bam_map_file(fh):
    """ Parse mapping of sample names to BAM files into a dict """
    return {name: bam for name, bam in map(lambda line: line.rstrip().split("\t"), fh)}


def parse_genes(genes_fh):
    """Return list of genes in genes file."""
    genes = []
    for line in genes_fh:
        genes.append(line.rstrip().split("\t")[0])
    return genes


def parse_vcf_file(vcf_reader, vep_cols, gene_list):
    """ Parse VCF file into a dict of records grouped by gene """
    records = defaultdict(list)
    for record in vcf_reader:
        vep_effect = parse_vep(vep_cols, record, tag="TOP_CSQ")[0]
        # If gene list not empty, skip genes not in list
        if gene_list and not (vep_effect["SYMBOL"] in gene_list or vep_effect["Gene"] in gene_list):
            continue
        gene = vep_effect["SYMBOL"] if vep_effect["SYMBOL"] else vep_effect["Gene"]
        records[gene].append(record)
    return records


def generate_autoigv_cmd(python_path, autoigv_path, genome):
    tmpl = "{} {} --file {} --directory {} --genome {} --mode 1 --host localhost --port 60151 --prefsfile ../{}\n"
    return tmpl.format(python_path, autoigv_path, POSITIONS_NAME, ".", genome, PREFS_NAME)


def generate_autoigv_positions(gene, records, bam_map):
    # Iterate over records
    lines = ""
    for record in records:
        lines += "{}:{}".format(record.CHROM, record.POS-1)
        for call in record.samples:
            if not call.is_variant:
                continue
            name = call.sample
            bam = bam_map[name]
            lines += "\t{}".format(bam)
        lines += "\n"
    return lines


def generate_master_script():
    script = "for GENE_DIR in $(find . -mindepth 1 -maxdepth 1 -type d)\ndo\n\t(cd $GENE_DIR && sh {})\ndone\n".format(SCRIPT_NAME)
    return script


def generate_prefs():
    prefs = "Order: host, port number, genome, default directory\nlocalhost\n60151\nhg19\n.\n"
    return prefs


if __name__ == '__main__':
    main()
