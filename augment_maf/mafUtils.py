# Useful functions for dealing with MAF files.

def get_protein_change(maf_row):
    """Parse the protein change from a MAF row"""
    try:
        prot_change = maf_row["Protein_Change"]
    except KeyError:
        prot_change = maf_row["HGVSp_Short"]
    return prot_change.replace("p.", "").replace("=", "")

def get_nref_allele(maf_row):
    """Get the nonreference allele from a MAF row"""
    ref_allele = maf_row["Reference_Allele"]
    if maf_row["Tumor_Seq_Allele1"] != ref_allele:
        return maf_row["Tumor_Seq_Allele1"]
    return maf_row["Tumor_Seq_Allele2"]
    
def is_snv(maf_row):
    """True if this row describes a SNV, False for an indel."""
    ref_allele = maf_row["Reference_Allele"]
    nref_allele = get_nref_allele(maf_row)
    return (len(ref_allele) == len(nref_allele) and 
            "-" not in [ref_allele, nref_allele])

def get_ensembl_id(maf_row):
    """Get the Ensembl gene ID for a MAF row."""
    if row["Entrez_Gene_Id"].startswith("ENSG"):
        return row["Entrez_Gene_Id"]
    return row["Gene"]

def get_base_change(maf_row):
    """Get the base change from a MAF row describing a SNV"""
    ref_allele = maf_row["Reference_Allele"]
    nref_allele = get_nref_allele(maf_row)
    return "{}>{}".format(ref_allele, nref_allele)

def get_cdna_change(maf_row):
    """Get the cDNA change from a MAF row describing a SNV"""
    try:
        return row["cDNA_Change"]
    except KeyError:
        cdna_ptn = "c[.](\d+)([ATCG])>([ATCG])"
        match = re.search(cdna_ptn, maf_row["HGVSc"])
        if match:
            return "{}{}{}".format(*match.group(2,1,3))
        return ""

def get_triplet(maf_row):
    """Get the triplet context from a MAF row describing a SNV"""
    try:
        return maf_row["Codons"].split("/")[0]
    except KeyError:
        return ""
    
def get_transcript(maf_row):
    """Get the transcript from a MAF row"""
    try:
        return row["Annotation_Transcript"]
    except KeyError:
        return row["Transcript_ID"]

def get_identifiers(maf_row):
    """Get all known identifiers from a MAF row"""
    try:
        ids = row["Existing_variation"].replace(",", "&")
    except KeyError:
        ids = row["dbSNP_RS"]
    ids = ids.replace(",", " ").replace(";", " ")
    ids = ids.split()
    return "&".join(id for id in ids if id != "novel")

def get_codon_pos(maf_row):
    """Get the codon position of a MAF row representing an indel"""
    try:
        positions = re.findall("\d+", maf_row["Transcript_Position"])
        positions = [int(int(i)/3) for i in positions]
    except KeyError:
        positions = [int(i) for i in re.findall("\d+", maf_row["HGVSc"])]

    if len(positions) == 1:
        return positions[0]
    else:
        return "{}-{}".format(min(positions), max(positions))
        
def get_effect(maf_row):
    """Get the effect from a MAF row representing an indel"""
    try:
        return maf_row["Consequence"]
    except KeyError:
        indel_class = maf_row["Variant_Classification"]
        if indel_class.startswith("Frame_Shift"):
            return "frameshift_variant"
        elif indel_class == "In_Frame_Del":
            return "inframe_deletion"
        elif indel_class == "In_Frame_Ins":
            return "inframe_insertion"
        else:
            msg = "Inserting NULL effect for an indel of type {}"
            warnings.warn(msg.format(indel_class))
            return ""

def get_allele_counts(maf_row):
    """Get allele counts from a row of a MAF file

    Returns a 4-tuple (normal_ref, normal_nref, tumour_ref, tumour_nref)"""
    counts = []
    for key in ["n_ref_count", "n_alt_count", "t_ref_count", "t_alt_count"]:
        try:
            counts.append(int(maf_row[key]))
        except (KeyError, ValueError):
            counts.append(0)
    return counts
