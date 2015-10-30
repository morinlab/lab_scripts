#!/usr/bin/env python

"""
annotate_vcf_on_cohort.py
=========================
This script is meant to filter several VCF files from
a given cohort simultaneously. This allows filtering
on the cohort level as opposed to simply on the sample
level. Effect prioritization based on vcf2maf.

Features
--------
- Filter out likely polymorphisms (seen in multiple
  individuals)

Dependencies
------------
- pyvcf

Known Issues
------------
- Assumes that all VCF files were annotated by VEP
  in the same way. This assumption is used for parsing
  the VEP columns from the header (taken from first file).
  It's also used when taking the consequence for a given
  variant from one of the samples.
- Assumes one sample per VCF file (first one)
"""

import sys
import argparse
import vcf
from vcf import utils
from vcf import model
import copy
import re
from collections import defaultdict

__version__ = "1.0.0"


def main():
    """Consolidate and annotate variants across cohort."""

    # Argument parsing
    args = parse_args()

    # Setup
    vcf_readers, vcf_iter = create_vcf_iter(*args.vcf)
    one_reader = copy.copy(vcf_readers[0])
    vep_cols = parse_vep_cols(one_reader)
    sample_names = [extract_sample_name(reader) for reader in vcf_readers]
    one_reader.samples = sample_names

    # Create a VCF Writer object with new INFO metadata lines
    one_reader_tweaked = tweak_vcf_reader(one_reader, args.snp_threshold, args.alt_codon_threshold, args.cluster_threshold)
    vcf_writer = vcf.Writer(args.output, template=one_reader_tweaked, lineterminator='\n')

    # Create dictionary of positions to annotated_snp_pos (empty if not specified)
    annotated_snp_pos_dict = create_pos_dict(args.annotated_snp_pos)

    # Create dictionary of positions in COSMIC
    cosmic_pos_dict = create_pos_dict(args.cosmic_pos)

    # Find mutated codons
    mutated_codons = find_mutated_codons(vcf_iter, vep_cols)

    # Detect potential hotspots in mutated_codons
    hotspot_codons = detect_hotspots(mutated_codons, args.alt_codon_threshold, args.cluster_threshold, args.max_for_hotspot)

    # Regenerate vcf_iter and annotate variants and output
    vcf_readers, vcf_iter = create_vcf_iter(*args.vcf)
    variant_iterator = annotate_variants(vcf_iter, vep_cols, sample_names, annotated_snp_pos_dict, cosmic_pos_dict, args.snp_threshold, hotspot_codons)
    for record in variant_iterator:
        vcf_writer.write_record(record)


def parse_args():
    """Parse the command-line arguments"""
    # Setup
    parser = argparse.ArgumentParser()
    # Optional arguments
    parser.add_argument("--output", "-o", default=sys.stdout, type=argparse.FileType("w"), help="Output VCF file")
    parser.add_argument("--snp_threshold", "-s", default=3, help="Min. number of cases for SNV -> SNP")
    parser.add_argument("--alt_codon_threshold", "-a", default=2, help="Min. number of alt. codons for being flagged as hotspot")
    parser.add_argument("--cluster_threshold", "-c", default=1, help="Max. distance between SNVs of a hotspot cluster")
    parser.add_argument("--annotated_snp_pos", type=argparse.FileType("r"), help="Annotated SNP positions (format: CHROM\\tPOS)")
    parser.add_argument("--max_for_hotspot", type=int, default=10, help="Max. number of cases for hotspot to be considered")
    parser.add_argument("--cosmic_pos", type=argparse.FileType("r"), help="COSMIC database of variants for annotation (format: CHROM\\tPOS)")
    # Positional arguments
    parser.add_argument("vcf", nargs="+", metavar="vcf_file")
    # Parsing
    args = parser.parse_args()
    return args


def create_vcf_iter(*vcf_paths):
    """Create independent VCF reader iterators.
    Returns a list of readers and a VCF iterator.
    """
    vcf_readers = [vcf.Reader(open(path)) for path in vcf_paths]
    vcf_iter = utils.walk_together(*vcf_readers, vcf_record_sort_key=lambda r: (r.CHROM, r.POS, r.REF, r.ALT))
    return vcf_readers, vcf_iter


def extract_sample_name(vcf_reader):
    """Extract name of the sample from a VCF reader.
    Only returns the first one.
    """
    vcf_record = vcf_reader.next()
    sample_name = vcf_record.samples[0].sample
    return sample_name


def parse_vep_cols(vcf_reader):
    """Parse VEP columns"""
    vep_desc = vcf_reader.infos["CSQ"].desc
    match = re.match(r".*Format: (.*)", vep_desc)
    vep_cols = match.group(1).split("|")
    return vep_cols


def parse_vep(vep_cols, vcf_record, tag="CSQ"):
    """Parse VEP INFO output and return as list"""
    vep_effects = []
    for vep_effect in vcf_record.INFO[tag]:
        vep_effect_vals = vep_effect.split("|")
        vep_effects.append(dict(zip(vep_cols, vep_effect_vals)))
    return vep_effects


def tweak_vcf_reader(vcf_reader, snp_threshold, alt_codon_threshold, cluster_threshold):
    """Tweak VCF reader to become VCF writer template"""
    for info in vcf_reader.infos:
        if info not in ["CSQ"]:
            del vcf_reader.infos[info]
    vcf_reader.infos["NUM_SAMPLES"] = vcf.parser._Info("NUM_SAMPLES", 1, "Integer", "Number of affected samples", __name__, __version__)
    vcf_reader.infos["TOP_CSQ"] = vcf.parser._Info("TOP_CSQ", ".", "String", "Top VEP effect", __name__, __version__)
    vcf_reader.infos["PROTEIN_CHANGE"] = vcf.parser._Info("PROTEIN_CHANGE", 0, "Flag", "Top effect changes protein", __name__, __version__)
    vcf_reader.infos["SNP"] = vcf.parser._Info("SNP", 0, "Flag", "Recurrent across population, namely in {} or more cases".format(snp_threshold), __name__, __version__)
    vcf_reader.infos["ANN_SNP_POS"] = vcf.parser._Info("ANN_SNP_POS", 0, "Flag", "Annotated SNP position", __name__, __version__)
    vcf_reader.infos["COSMIC"] = vcf.parser._Info("COSMIC", 0, "Flag", "Annotated COSMIC position", __name__, __version__)
    vcf_reader.infos["HOTSPOT"] = vcf.parser._Info("HOTSPOT", 0, "Flag", "Codon altered in {}+ ways".format(alt_codon_threshold), __name__, __version__)
    vcf_reader.infos["HOTSPOT_CLUSTER"] = vcf.parser._Info("HOTSPOT_CLUSTER", 0, "Flag", "Codons within {} positions are mutated".format(cluster_threshold), __name__, __version__)
    return vcf_reader


def create_pos_id(chrom, pos):
    """Create unique ID for chromosome and position"""
    identifier = "{}_{}".format(chrom, pos)
    return identifier


def create_pos_dict(pos_file):
    """Create dictionary of positions for fast lookup."""
    positions = {}
    if pos_file is None:
        return positions
    for line in pos_file:
        if line.startswith("#"):
            continue
        chrom, pos = line.split("\t")[0:2]
        identifier = create_pos_id(chrom, pos)
        positions[identifier] = 1
    return positions


def prioritize_effects(vep_effects):
    """Prioritize VEP effects. Modelled after vcf2maf implementation:
    https://github.com/mskcc/vcf2maf
    But with differences. Priotized as follows:
    1) Transcript biotype priority
    2) Consequence severity
    3) Canonical transcript?
    4) Decreasing transcript length
    """
    all_priorities = []  # List of tuples corresponding to the three priorities
    for effect in vep_effects:
        # Transcript biotype priority
        biotype = effect["BIOTYPE"]
        biotype_priority = BIOTYPE_PRIORITY[biotype]
        # Consequence severity priority
        consequence = effect["Consequence"]
        if "&" in consequence:
            consequences = consequence.split("&")
            consequence_priority = min([CONSEQUENCE_PRIORITY[c] for c in consequences])
        else:
            consequence_priority = CONSEQUENCE_PRIORITY[consequence]
        # Canonical transcript prioritization
        # Since higher priority is a lower number, True -> 1 and False -> 2
        if effect["CANONICAL"] == "YES":
            canonical_priority = 1
        else:
            canonical_priority = 2
        # Transcript length priority
        # Since length is inversely proportional to priotity, convert to negative value
        cdna_info = effect["cDNA_position"]
        if cdna_info == "":
            transcript_length = 0
        else:
            match = re.search(r'\/(\d+)$', cdna_info)
            transcript_length = int(match.group(1))
        transcript_length_priority = 0 - transcript_length
        # Combine priorities
        combined = (biotype_priority, consequence_priority, canonical_priority, transcript_length_priority)
        all_priorities.append(combined)
    combined = zip(vep_effects, all_priorities)
    combined_sorted = sorted(combined, key=lambda x: x[1])
    ranked_effects = [eff[0] for eff in combined_sorted]
    return ranked_effects


def obtain_one_record(records):
    """Obtain one exemplary record"""
    one_record = None
    for record in records:
        if record is None:
            continue
        one_record = record
        break
    return one_record


def obtain_top_effect(vcf_record, vep_cols):
    """Obtain top VEP effect"""
    effects = parse_vep(vep_cols, vcf_record)
    ranked_effects = prioritize_effects(effects)
    top_effect = ranked_effects[0]
    return top_effect


def obtain_mutated_codon_info(vep_effect):
    """Obtain information about mutated codon.
    Returns tuple: (transcript_name, codon_position, alternate_codon)
    """
    transcript = vep_effect["Feature"]
    codon_num = int(re.match(r"(\d+)", vep_effect["Protein_position"]).group(1))
    alt_codon = re.match(r"\w{3}/(\w{3})", vep_effect["Codons"]).group(1)
    return transcript, codon_num, alt_codon


def find_mutated_codons(vcf_iter, vep_cols):
    """Find mutated codons and returns a dictionary
    of transcripts as keys and dictionaries as values.
    These dictionaries in turn contain codon numbers as
    keys and dicts of alternate codon sequences as values.
    These dicts contain the number of instances this
    alteration is seen in the cohort.
    See below for example.
        {'ENSCAFT00000001435': {
            64: {'tAc': 2},
            71: {'Acc': 4} }}
    To be used in tandem with detect_hotspots().
    """
    mutated_codons = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for records in vcf_iter:
        # Find one exemplary record/sample carrying variant
        one_record = obtain_one_record(records)
        # Obtain top effect
        top_effect = obtain_top_effect(one_record, vep_cols)
        # If SNV and non-synonymous, increment variant codon in mutated_codons by num_samples
        if one_record.is_snp and top_effect["HGVSp"] != "":
            transcript, codon_num, alt_codon = obtain_mutated_codon_info(top_effect)
            num_samples = sum([1 for r in records if r is not None])
            mutated_codons[transcript][codon_num][alt_codon] += num_samples
    return mutated_codons


def detect_hotspots(mutated_codons, alt_codon_threshold, cluster_threshold, max_for_hotspot):
    """Detects hotspot codons when given a dictionary
    of transcripts as keys and dictionaries as values.
    These dictionaries in turn contain codon numbers as
    keys and sets of alternate codon sequences as values.
    See below for example.
        {'ENSCAFT00000001435': {
            64: {'tAc': 2},
            71: {'Acc': 4} }}
    Returns a dictionary of trancripts-dictionary key-value
    pairs, where the dictionary contains the codon numbers
    as keys and whether it is flagged as a hotspot as values.
    See below for example.
        {'ENSCAFT00000001435': {
            64: True,
            71: False }}
    """
    hotspot_codons = {}
    for transcript, codons_dict in mutated_codons.items():
        # Initialize hotspot_codons
        hotspot_codons[transcript] = {codon_num: False for codon_num in codons_dict.keys()}
        # Identify codons altered in different ways
        for codon_num, alt_codons_dict in codons_dict.items():
            if len(alt_codons_dict) >= alt_codon_threshold and not any(map(lambda num: num > max_for_hotspot, alt_codons_dict.values())):
                hotspot_codons[transcript][codon_num] = True
        # Determine clustered mutated codons
        # ordered_codons = sorted(codons_dict.keys())
        # prev_codon = None
        # for curr_codon in ordered_codons:
        #     if prev_codon is not None and (curr_codon - prev_codon <= cluster_threshold):
        #         hotspot_codons[transcript][curr_codon] = True
        #         hotspot_codons[transcript][prev_codon] = True
        #     prev_codon = curr_codon
    return hotspot_codons


def annotate_variants(vcf_iter, vep_cols, sample_names, annotated_snp_pos_dict, cosmic_pos_dict, snp_threshold, hotspot_codons):
    """Generator yielding annotated records ready
    for being outputted to a file. Combine samples
    together.
    """
    for records in vcf_iter:
        # Find one exemplary record/sample carrying variant
        one_record = obtain_one_record(records)
        # Otherwise, continue
        # Find common denominator of format fields
        all_fmt_fields = [set(r.FORMAT.split(":")) for r in records if r is not None]
        common_fmt_fields = set.intersection(*all_fmt_fields)
        # Ensure that GT is first
        if "GT" in common_fmt_fields:
            common_fmt_fields = ["GT"] + list(common_fmt_fields - set(["GT"]))
        common_fmt_string = ":".join(common_fmt_fields)
        one_record.FORMAT = common_fmt_string
        format_tuple = model.make_calldata_tuple(common_fmt_fields)
        blanks = ["0"] * (len(common_fmt_fields) - 1)
        uncalled_format = format_tuple("0/0", *blanks)
        # Obtain top effect
        top_effect = obtain_top_effect(one_record, vep_cols)
        # Generate calls for all records
        combined_samples = []
        num_affected_samples = 0
        for sample_name, record in zip(sample_names, records):
            if record is None:
                new_sample = vcf.model._Call(site=one_record, sample=sample_name, data=uncalled_format)
            else:
                new_sample = record.samples[0]
                fmt_dict = record.samples[0].data._asdict()
                fmt_data = format_tuple(*[fmt_dict[field] for field in common_fmt_fields])
                new_sample.data = fmt_data
                num_affected_samples += 1
            combined_samples.append(new_sample)
        one_record.samples = combined_samples
        # Remove unnecessary INFO tags
        for k in one_record.INFO.keys():
            if k not in []:
                del one_record.INFO[k]
        # Clear QUAL col
        one_record.QUAL = "."
        # Annotate variant
        variant_id = create_pos_id(one_record.CHROM, one_record.POS)
        if variant_id in annotated_snp_pos_dict:
            one_record.INFO["ANN_SNP_POS"] = True
        if variant_id in cosmic_pos_dict:
            one_record.INFO["COSMIC"] = True
        one_record.INFO["NUM_SAMPLES"] = num_affected_samples
        top_effect_encoded = "|".join([top_effect[col] for col in vep_cols])
        one_record.INFO["TOP_CSQ"] = top_effect_encoded
        if num_affected_samples >= snp_threshold:
            one_record.INFO["SNP"] = True
        if top_effect["HGVSp"] != "" and "synonymous_variant" not in top_effect["Consequence"]:
            one_record.INFO["PROTEIN_CHANGE"] = True
            if one_record.is_snp:
                transcript, codon_num, alt_codon = obtain_mutated_codon_info(top_effect)
                if hotspot_codons[transcript].setdefault(codon_num, False):
                    one_record.INFO["HOTSPOT"] = True
        # Yield variant
        yield one_record


BIOTYPE_PRIORITY = {
    'protein_coding': 1,  # Contains an open reading frame (ORF)
    'LRG_gene': 2,  # Gene in a "Locus Reference Genomic" region known to have disease-related sequence variations
    'IG_C_gene': 2,  # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
    'IG_D_gene': 2,  # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
    'IG_J_gene': 2,  # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
    'IG_V_gene': 2,  # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
    'TR_C_gene': 2,  # T-cell receptor (TcR) genes imported or annotated according to the IMGT
    'TR_D_gene': 2,  # T-cell receptor (TcR) genes imported or annotated according to the IMGT
    'TR_J_gene': 2,  # T-cell receptor (TcR) genes imported or annotated according to the IMGT
    'TR_V_gene': 2,  # T-cell receptor (TcR) genes imported or annotated according to the IMGT
    'miRNA': 3,  # Non-coding RNA predicted using sequences from RFAM and miRBase
    'snRNA': 3,  # Non-coding RNA predicted using sequences from RFAM and miRBase
    'snoRNA': 3,  # Non-coding RNA predicted using sequences from RFAM and miRBase
    'rRNA': 3,  # Non-coding RNA predicted using sequences from RFAM and miRBase
    'lincRNA': 3,  # Long, intervening noncoding (linc) RNAs, that can be found in evolutionarily conserved, intergenic regions
    'Mt_tRNA': 4,  # Non-coding RNA predicted using sequences from RFAM and miRBase
    'Mt_rRNA': 4,  # Non-coding RNA predicted using sequences from RFAM and miRBase
    'antisense': 5,  # Has transcripts that overlap the genomic span (i.e. exon or introns) of a protein-coding locus on the opposite strand
    'sense_intronic': 5,  # Long non-coding transcript in introns of a coding gene that does not overlap any exons
    'sense_overlapping': 5,  # Long non-coding transcript that contains a coding gene in its intron on the same strand
    '3prime_overlapping_ncrna': 5,  # Transcripts where ditag and/or published experimental data strongly supports the existence of short non-coding transcripts transcribed from the 3'UTR
    'misc_RNA': 5,  # Non-coding RNA predicted using sequences from RFAM and miRBase
    'non_coding': 5,  # Transcript which is known from the literature to not be protein coding
    'regulatory_region': 6,  # A region of sequence that is involved in the control of a biological process
    'disrupted_domain': 6,  # Otherwise viable coding region omitted from this alternatively spliced transcript because the splice variation affects a region coding for a protein domain
    'processed_transcript': 6,  # Doesn't contain an ORF
    'TEC': 6,  # To be Experimentally Confirmed. This is used for non-spliced EST clusters that have polyA features. This category has been specifically created for the ENCODE project to highlight regions that could indicate the presence of protein coding genes that require experimental validation, either by 5' RACE or RT-PCR to extend the transcripts, or by confirming expression of the putatively-encoded peptide with specific antibodies
    'TF_binding_site': 7,  # A region of a nucleotide molecule that binds a Transcription Factor or Transcription Factor complex
    'CTCF_binding_site': 7,  # A transcription factor binding site with consensus sequence CCGCGNGGNGGCAG, bound by CCCTF-binding factor
    'promoter_flanking_region': 7,  # A region immediately adjacent to a promoter which may or may not contain transcription factor binding sites
    'enhancer': 7,  # A cis-acting sequence that increases the utilization of (some) eukaryotic promoters, and can function in either orientation and in any location (upstream or downstream) relative to the promoter
    'promoter': 7,  # A regulatory_region composed of the TSS(s) and binding sites for TF_complexes of the basal transcription machinery
    'open_chromatin_region': 7,  # A DNA sequence that in the normal state of the chromosome corresponds to an unfolded, un-complexed stretch of double-stranded DNA
    'retained_intron': 7,  # Alternatively spliced transcript believed to contain intronic sequence relative to other, coding, variants
    'nonsense_mediated_decay': 7,  # If the coding sequence (following the appropriate reference) of a transcript finishes >50bp from a downstream splice site then it is tagged as NMD. If the variant does not cover the full reference coding sequence then it is annotated as NMD if NMD is unavoidable i.e. no matter what the exon structure of the missing portion is the transcript will be subject to NMD
    'non_stop_decay': 7,  # Transcripts that have polyA features (including signal) without a prior stop codon in the CDS, i.e. a non-genomic polyA tail attached directly to the CDS without 3' UTR. These transcripts are subject to degradation
    'ambiguous_orf': 7,  # Transcript believed to be protein coding, but with more than one possible open reading frame
    'pseudogene': 8,  # Have homology to proteins but generally suffer from a disrupted coding sequence and an active homologous gene can be found at another locus. Sometimes these entries have an intact coding sequence or an open but truncated ORF, in which case there is other evidence used (for example genomic polyA stretches at the 3' end) to classify them as a pseudogene. Can be further classified as one of the following
    'processed_pseudogene': 8,  # Pseudogene that lack introns and is thought to arise from reverse transcription of mRNA followed by reinsertion of DNA into the genome
    'polymorphic_pseudogene': 8,  # Pseudogene owing to a SNP/DIP but in other individuals/haplotypes/strains the gene is translated
    'retrotransposed': 8,  # Pseudogene owing to a reverse transcribed and re-inserted sequence
    'translated_processed_pseudogene': 8,  # Pseudogenes that have mass spec data suggesting that they are also translated
    'translated_unprocessed_pseudogene': 8,  # Pseudogenes that have mass spec data suggesting that they are also translated
    'transcribed_processed_pseudogene': 8,  # Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
    'transcribed_unprocessed_pseudogene': 8,  # Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
    'unitary_pseudogene': 8,  # A species specific unprocessed pseudogene without a parent gene, as it has an active orthologue in another species
    'unprocessed_pseudogene': 8,  # Pseudogene that can contain introns since produced by gene duplication
    'Mt_tRNA_pseudogene': 8,  # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    'tRNA_pseudogene': 8,  # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    'snoRNA_pseudogene': 8,  # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    'snRNA_pseudogene': 8,  # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    'scRNA_pseudogene': 8,  # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    'rRNA_pseudogene': 8,  # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    'misc_RNA_pseudogene': 8,  # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    'miRNA_pseudogene': 8,  # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    'IG_C_pseudogene': 8,  # Inactivated immunoglobulin gene
    'IG_J_pseudogene': 8,  # Inactivated immunoglobulin gene
    'IG_V_pseudogene': 8,  # Inactivated immunoglobulin gene
    'TR_J_pseudogene': 8,  # Inactivated immunoglobulin gene
    'TR_V_pseudogene': 8,  # Inactivated immunoglobulin gene
    'artifact': 9,  # Used to tag mistakes in the public databases (Ensembl/SwissProt/Trembl)
    '': 9
}


CONSEQUENCE_PRIORITY = {
    'transcript_ablation': 1,  # A feature ablation whereby the deleted region includes a transcript feature
    'exon_loss_variant': 1,  # A sequence variant whereby an exon is lost from the transcript
    'splice_donor_variant': 2,  # A splice variant that changes the 2 base region at the 5' end of an intron
    'splice_acceptor_variant': 2,  # A splice variant that changes the 2 base region at the 3' end of an intron
    'stop_gained': 3,  # A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript
    'frameshift_variant': 3,  # A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three
    'stop_lost': 3,  # A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript
    'start_lost': 4,  # A codon variant that changes at least one base of the canonical start codon
    'initiator_codon_variant': 4,  # A codon variant that changes at least one base of the first codon of a transcript
    'disruptive_inframe_insertion': 5,  # An inframe increase in cds length that inserts one or more codons into the coding sequence within an existing codon
    'disruptive_inframe_deletion': 5,  # An inframe decrease in cds length that deletes bases from the coding sequence starting within an existing codon
    'inframe_insertion': 5,  # An inframe non synonymous variant that inserts bases into the coding sequence
    'inframe_deletion': 5,  # An inframe non synonymous variant that deletes bases from the coding sequence
    'missense_variant': 6,  # A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved
    'conservative_missense_variant': 6,  # A sequence variant whereby at least one base of a codon is changed resulting in a codon that encodes for a different but similar amino acid. These variants may or may not be deleterious
    'rare_amino_acid_variant': 6,  # A sequence variant whereby at least one base of a codon encoding a rare amino acid is changed, resulting in a different encoded amino acid
    'transcript_amplification': 7,  # A feature amplification of a region containing a transcript
    'stop_retained_variant': 8,  # A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
    'synonymous_variant': 8,  # A sequence variant where there is no resulting change to the encoded amino acid
    'splice_region_variant': 9,  # A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron
    'incomplete_terminal_codon_variant': 10,  # A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed
    'protein_altering_variant': 11,  # A sequence variant which is predicted to change the protein encoded in the coding sequence
    'coding_sequence_variant': 11,  # A sequence variant that changes the coding sequence
    'mature_miRNA_variant': 11,  # A transcript variant located with the sequence of the mature miRNA
    'exon_variant': 11,  # A sequence variant that changes exon sequence
    '5_prime_UTR_variant': 12,  # A UTR variant of the 5' UTR
    '5_prime_UTR_premature_start_codon_gain_variant': 12,  # snpEff-specific effect, creating a start codon in 5' UTR
    '3_prime_UTR_variant': 12,  # A UTR variant of the 3' UTR
    'non_coding_exon_variant': 13,  # A sequence variant that changes non-coding exon sequence
    'non_coding_transcript_exon_variant': 13,  # snpEff-specific synonym for non_coding_exon_variant
    'non_coding_transcript_variant': 14,  # A transcript variant of a non coding RNA gene
    'nc_transcript_variant': 14,  # A transcript variant of a non coding RNA gene (older alias for non_coding_transcript_variant)
    'intron_variant': 14,  # A transcript variant occurring within an intron
    'intragenic_variant': 14,  # A variant that occurs within a gene but falls outside of all transcript features. This occurs when alternate transcripts of a gene do not share overlapping sequence
    'INTRAGENIC': 14,  # snpEff-specific synonym of intragenic_variant
    'NMD_transcript_variant': 15,  # A variant in a transcript that is the target of NMD
    'upstream_gene_variant': 16,  # A sequence variant located 5' of a gene
    'downstream_gene_variant': 16,  # A sequence variant located 3' of a gene
    'TFBS_ablation': 17,  # A feature ablation whereby the deleted region includes a transcription factor binding site
    'TFBS_amplification': 17,  # A feature amplification of a region containing a transcription factor binding site
    'TF_binding_site_variant': 17,  # A sequence variant located within a transcription factor binding site
    'regulatory_region_ablation': 17,  # A feature ablation whereby the deleted region includes a regulatory region
    'regulatory_region_amplification': 17,  # A feature amplification of a region containing a regulatory region
    'regulatory_region_variant': 17,  # A sequence variant located within a regulatory region
    'regulatory_region': 17,  # snpEff-specific effect that should really be regulatory_region_variant
    'feature_elongation': 18,  # A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence
    'feature_truncation': 18,  # A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence
    'intergenic_variant': 19,  # A sequence variant located in the intergenic region, between genes
    'intergenic_region': 19,  # snpEff-specific effect that should really be intergenic_variant
    '': 20
}


if __name__ == "__main__":
    main()
