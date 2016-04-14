import argparse
parser = argparse.ArgumentParser(description='parse bed files and search for recurrently affected genes and artifacts')


parser.add_argument('-b','--bed_files', action="store", dest="bed",nargs="+")
parser.add_argument('-g','--gene_list_file',action="store",dest="gene_file")
parser.add_argument('-c','--case_list_file',action="store",dest="case_file",help="file containing cases of interest")


args = parser.parse_args()

bed = args.bed
known_genes = []

gene_file = args.gene_file
if gene_file:
    gene_h = open(gene_file,"r")
    for line in gene_h:
        line = line.rstrip("\n")
        known_genes.append(line)


selected_cases = []
case_file = args.case_file
if case_file:
    case_h = open(case_file,"r")
    for line in case_h:
        line = line.rstrip("\n")
        selected_cases.append(line)



gene_groups = []
print "%s bed files" % len(bed)

def bed2string(bed_data):
    col_ids = ["CHROM_A","START_A","END_A","CHROM_B","START_B","END_B","ID","QUAL","STRAND_A","STRAND_B","TYPE","FILTER","INFO","FORMAT","detail1","detail2","genechrom","genestart","geneend","gene","file"]
    out_s = ""

    for col in col_ids:
        out_s = out_s + "%s\t" % bed_data[col]
    return(out_s)
def parseBedpeLine(line):
    cols = line.split()
    ##CHROM_A START_A END_A CHROM_B START_B END_B ID QUAL STRAND_A STRAND_B TYPE FILTER INFO FORMAT
    col_ids = ["CHROM_A","START_A","END_A","CHROM_B","START_B","END_B","ID","QUAL","STRAND_A","STRAND_B","TYPE","FILTER","INFO","FORMAT","detail1","detail2","genechrom","genestart","geneend","gene"]
    parsed = {}
    i = 0
    for col in cols:
        try:
            parsed[col_ids[i]] = col
        except IndexError:
            pass
        i+=1
    parsed["START_A"] = int(parsed["START_A"])
    return parsed

def updateGeneGroups(genes):
    '''Take a list of genes and either merge with group containing one or more shared genes or create a new group'''
    found = 0
    for gene in genes:
        group_index = 0
        for group in gene_groups:
            if gene in group:
                for gene in genes:
                    group.add(gene)
                return group_index
            group_index+=1
    #no group for any gene, start new group
    new_group = set(genes)

    gene_groups.append(new_group)

    return(len(gene_groups)-1)

gene_detail = {} #key on individual genes (eventually switch to gene groups)
gene_count = {}

event_selected_cases = {} #for each event/gene group, track (using a set) all the "selected" cases that the event was seen in
event_unselected_cases = {}
last_start = 0
this_gene_list = []
this_group = 0
saved_bedcols = []
bed_count = 0
for bed_file in bed:
    bed_h = open(bed_file,"r")
    bed_count+=1
    gene_good = {}
    gene_bad = {}
    group_counted = {}
    print bed_file
    for line in bed_h:
        line = line.rstrip("\n")
        current_bedcols = parseBedpeLine(line)
        current_bedcols["file"] = bed_file
        
        if last_start == current_bedcols["START_A"]:
            #same breakpoint, adding another gene
            try:
                (gene,ensg) = current_bedcols['gene'].split("_")
            except ValueError:
                gene = current_bedcols['gene']
            this_gene_list.append(gene)
            saved_bedcols.append(current_bedcols)
        else:
            #update groups with the genes in this list
            if len(this_gene_list):
                #print "UPDATING WITH:" 
                #print this_gene_list
                this_group = updateGeneGroups(this_gene_list)
                if len(selected_cases):
                    if len(saved_bedcols):
                        filename = saved_bedcols[0]['file']
                    sel = 0
                    for case in selected_cases:
                        if case in filename:
                            sel = 1
                            if event_selected_cases.has_key(this_group):
                                event_selected_cases[this_group].add(case)
                            else:
                                event_selected_cases[this_group] = set([case])
                    if not sel:
                        if event_unselected_cases.has_key(this_group):
                            event_unselected_cases[this_group].add(filename)
                        else:
                            event_unselected_cases[this_group] = set([filename])

                passed = 0
                for bedcols in saved_bedcols:
                    if bedcols["FILTER"] == "PASS":
                        passed = 1
                    if gene_detail.has_key(this_group):
                        gene_detail[this_group].append(bedcols)
                        
                    else:
                        gene_detail[this_group] = [bedcols]
                        #gene_count[this_group] = {bedcols["FILTER"]:1}
                    break #only count one of the duplicate bed lines
                if not group_counted.has_key(this_group):
                    group_counted[this_group] = 1
                    if passed:
                        gene_good[this_group] = 1
                    else:
                        gene_bad[this_group] = 1
                    if gene_count.has_key(this_group):
                        if passed:
                            try:
                                gene_count[this_group]["PASS"]+=1
                            except KeyError:
                                gene_count[this_group]["PASS"] = 1
                        else:
                            try:
                                gene_count[this_group]["MinSomaticScore"]+=1
                            except KeyError:
                                gene_count[this_group]["MinSomaticScore"] = 1
                    else:
                        if passed:
                            gene_count[this_group]  = {"PASS":1}
                        else:
                            gene_count[this_group] = {"MinSomaticScore":1}
            try:
                (gene,ensg) = current_bedcols['gene'].split("_")
            except ValueError:
                gene = current_bedcols['gene']
            this_gene_list = [gene]
            last_start = current_bedcols["START_A"]
            saved_bedcols = []

out_h = open("gene_group_summary.txt","w")
for group in gene_detail.keys():
    for bed_data in gene_detail[group]:
        string = bed2string(bed_data)
        string = str(group) + "\t" + string
        out_h.write(string)
        relevant_genes = []
        for gene in gene_groups[group]:
            if gene in known_genes:
                relevant_genes.append(gene)
        if len(relevant_genes):
            #print "YAY"
            #print relevant_genes
            #exit()
            out_h.write(str(relevant_genes))
        else:
            out_h.write(str(gene_groups[group]))
        out_h.write("\n")
out_h.close()


for gene in gene_count.keys():
    good = 0
    filt = 0
    try:
        good = gene_count[gene]["PASS"]
    except KeyError:
        pass
    try:
        filt = gene_count[gene]["MinSomaticScore"]
    except KeyError:
        pass
    relevant_genes = []
    for genes in gene_groups[gene]:
        if genes in known_genes:
            relevant_genes.append(genes)
    selected_freq = 0
    unselected_freq = 0
    if event_selected_cases.has_key(gene):
        selected_freq = len(event_selected_cases[gene])
    if event_unselected_cases.has_key(gene):
        unselected_freq =len(event_unselected_cases[gene])
    num_sel = len(selected_cases)
    num_unsel = bed_count - num_sel
    perc_sel = 100* float(selected_freq) /float(num_sel)
    perc_unsel = 100 * float(unselected_freq) / float(num_unsel)

    if len(relevant_genes):
        if good + filt > 3:
            good_ratio = float(good+1)/float(filt+1)
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (gene,good,filt,selected_freq,perc_sel,unselected_freq,perc_unsel,str(relevant_genes))
    else:
        if good + filt > 3:
            good_ratio = float(good+1)/float(filt+1)
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (gene,good,filt,selected_freq,perc_sel,unselected_freq,perc_unsel,gene_groups[gene])
