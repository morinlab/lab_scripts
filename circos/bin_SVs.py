#!/usr/bin/env python
import vcf
import math
import glob, os


fai = "/projects/rmorin/common/genomes/human_all.fasta.fai" #for hg18
#fai = '/projects/rmorin/common/genomes/hg19/GRCh37-lite.fa.fai'
#vcf_path = '/projects/rmorin/projects/nhl_meta_analysis/results/Manta/DLBCL_refractory_Gascoyne/'
vcf_path = '/projects/rmorin/projects/nhl_meta_analysis/results/Manta/DLBCL_Morin/'
chroms = range(1,23)
chroms.append("X")

score_cutoff = 25

#chroms = (14,18)




def main(fai_file,vcf_file_location,bin_size=5000000):
    bins = createGenomeBins(fai,bin_size)
    bin_break_counts = {} #store a count for every bin
    bin_break_pairs = {} #store the pair details for every event, organized by bin for retrieval
    outh = open("summary_svs.txt","w")
    os.chdir(vcf_file_location)
    sample = ""
    for file in glob.glob("*.vcf"):
        sample = file.rstrip(".vcf")
        print "%s for %s"% (sample,file)
        file = vcf_path + file
        print(file)
        tallyBins(file,bins,bin_break_counts,bin_break_pairs,outh,sample)
    print bin_break_pairs.keys()
    #print bin_break_counts.keys()
    #exit()
    #print bin_break_pairs
    for chrom in bin_break_pairs.keys():
        hit_bins = bin_break_counts[chrom].keys()
        for hit_bin in hit_bins:

            count = bin_break_counts[chrom][hit_bin]
            if count <= 2:
                continue
            start = bins[chrom][0][hit_bin]
            end = bins[chrom][1][hit_bin]

            if bin_break_pairs.has_key(chrom):
                for bin in bin_break_pairs[chrom].keys():
                    if bin == hit_bin:
                        dup_count = {}
                        #seen_before = {}
                        for linkage in bin_break_pairs[chrom][bin]:
                            (c1,c2,b1,b2,sample,type) = linkage
                            #if c1 == c2:
                            #    print "INTRACHROMOSOMAL"
                            #    print linkage
                            #    exit()
                            if dup_count.has_key(sample):
                                continue
                            dup_count[file] = 1
                            
                            s1 = bins[c1][0][b1]
                            e1 = bins[c1][1][b1]
                            s2 = bins[c2][0][b2]
                            e2 = bins[c2][1][b2]
                            print "%s\tchr%s\tchr%s\t%s\t%s\t%s" % (sample,c1,c2,type,s1,s2) 
            else:
                print "no chromosome for %s" % chrom
    outh.close()
def getBin(starts,ends,pos):
    for i in range(0,len(ends)):
        start = starts[i]
        end = ends[i]
        if pos >= start and pos <= end:
            return(i)

    print "%s NOT FOUND!!!" % pos
    print ends[-1]
    exit()

def tallyBins(vcf_file,bins,tally,links,outh,sample):
    #first load chrom, pos details for all and key on ID
    sv_id = {}
    sv_mates = {}
    vcf_reader = vcf.Reader(open(vcf_file,"r"))
    for record in vcf_reader:
        sv_id[record.ID] = record
        if record.INFO.has_key("MATEID"):
            mate_id = record.INFO["MATEID"][0]
            sv_mates[record.ID] = mate_id

    for sv in sv_mates.keys():
        chrom1 = sv_id[sv].CHROM
        chrom2 = sv_id[sv_mates[sv]].CHROM
        #if chrom1 == chrom2:
        #    print "INTRACHROMOSOMAL"
        #    print sv_id[sv]
        #    print sv_id[sv_mates[sv]]
    
    for chrom in bins.keys():
        #print "===============%s<<<<<<<<<<<<<<<<<<" % chrom
        
        vcf_reader = vcf.Reader(open(vcf_file,"r"))
        for record in vcf_reader:
            if not record.CHROM == chrom:
                continue
            if record.INFO["SVTYPE"] == "INV":
                #both breaks are on the same chromosome, find bin for each
                inv_end = record.INFO["END"]
                score = record.INFO["SOMATICSCORE"]
                if score < score_cutoff:
                    continue
                string = "%s %s %s INV %s\n" % (chrom,record.POS,inv_end,vcf_file) 
                outh.write(string)
                bin_inds = getBin(bins[chrom][0],bins[chrom][1],record.POS)
                bin_inde = getBin(bins[chrom][0],bins[chrom][1],inv_end)
                #print "binS: %s binE: %s for inversion %s %s" % (bin_inds,bin_inde,record.POS,inv_end)
                if tally.has_key(chrom):
                    if tally[chrom].has_key(bin_inds):
                        tally[chrom][bin_inds] +=1

                    else:
                        tally[chrom][bin_inds] = 1
                else:
                    tally[chrom] = {bin_inds:1}
                
                if tally.has_key(chrom):
                    if tally[chrom].has_key(bin_inde):
                        tally[chrom][bin_inde] +=1

                    else:
                        tally[chrom][bin_inde] = 1
                else:
                    tally[chrom] = {bin_inde:1}

                link = (chrom,chrom,bin_inds,bin_inde,sample,"INV")
                if links.has_key(link[0]):
                    if links[link[0]].has_key(link[2]):
                        links[link[0]][link[2]].append(link)
                    else:
                        links[link[0]][link[2]] = [link]
                else:
                    links[link[0]] = {link[2]:[link]}
            #continue
            if record.INFO.has_key("MATEID"):
                mate_id = record.INFO["MATEID"][0]
                score = record.INFO["SOMATICSCORE"]
                if score < score_cutoff:
                    continue

                mate_record = sv_id[mate_id]
                bin_ind = getBin(bins[chrom][0],bins[chrom][1],record.POS)
                mate_chr = mate_record.CHROM
                ints = 0
                try:
                    ints = int(mate_chr)
                except ValueError:
                    pass

                if mate_chr in chroms or ints in chroms:
                    pass
                else:
                    #print "lacking %s, %s" % (mate_chr,ints)
                    continue
                mate_bin_ind = getBin(bins[mate_chr][0],bins[mate_chr][1],mate_record.POS)
                twochrom = [mate_chr,chrom]
                twochrom.sort()
                link = twochrom
                #print twochrom
                if twochrom[0] == mate_chr:
                    link.append(mate_bin_ind)
                    link.append(bin_ind)
                else:
                    link.append(bin_ind)
                    link.append(mate_bin_ind)
                link.append(sample)
                link.append("TRA")
                string = "%s %s %s %s TRA %s\n" % (chrom,mate_chr,record.POS,mate_record.POS,vcf_file)
                outh.write(string)

                #print "using %s" % link[0]
                if links.has_key(link[0]):
                    if links[link[0]].has_key(link[2]):
                        links[link[0]][link[2]].append(link)
                    else:
                        links[link[0]][link[2]] = [link]
                else:
                    links[link[0]] = {link[2]:[link]}
                #if twochrom[0] == "X" or twochrom[1] == "X":
                #    print links[link[0]]
                #    exit()
                #print "This: %s %s %s (%s) Mate: %s %s %s (%s)" % (chrom,record.POS,bin_ind,bins[chrom][0][bin_ind],mate_chr,mate_record.POS,mate_bin_ind,bins[mate_chr][0][mate_bin_ind])
                #exit()
                #print "%s in bin %s: %s-%s" % (record.POS,bin_ind,bins[chrom][0][bin_ind],bins[chrom][1][bin_ind])
                if tally.has_key(chrom):
                    if tally[chrom].has_key(bin_ind):
                        tally[chrom][bin_ind] +=1
                        
                    else:
                        tally[chrom][bin_ind] = 1
                else:
                    tally[chrom] = {bin_ind:1}

        #(bin_starts,bin_ends) = bins[chrom]
        #for i in range(0,len(bin_ends)):                                                                                                                  #     for record in vcf_reader:   
        #        print "%s %s %s-%s" % (i,vals[0],bin_starts[i],bin_ends[i])
        #        print record
        

def createGenomeBins(fai,bin_size):
    chrom_bins = {}
    fai_handle = open(fai,"r")
    for line in fai_handle:
        line = line.rstrip("\n")
        vals=line.split("\t")
        ints = 0
        try:
            ints = int(vals[0])
        except ValueError:
            pass
        if vals[0] in chroms or ints in chroms:
            #break chromosome into bins 
            chromlen = float(vals[1])
            num_bins = int(math.floor(chromlen / bin_size))
            #print "need %s bins" % num_bins
            bin_starts = range(1,int(chromlen),bin_size)
            bin_ends = bin_starts[1:len(bin_starts)]
            bin_ends.append(bin_ends[-1]+bin_size)
            bin_ends = map(lambda x: x-1, bin_ends)
            #for i in range(0,len(bin_ends)):
            #    print "%s %s %s-%s" % (i,vals[0],bin_starts[i],bin_ends[i])
            chrom_bins[vals[0]] = (bin_starts,bin_ends)
    return(chrom_bins)

main(fai,vcf_path)
