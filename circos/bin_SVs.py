#!/usr/bin/env python
import vcf
import math
import glob, os
import argparse

chroms = range(1,23)
chroms.append("X")

score_cutoff = 25

#chroms = (14,18)




def main():
    parser = argparse.ArgumentParser(description="summarize intra- and intrachromosomal SV information from a series of vcf files")
    parser.add_argument("--fai", "-f",dest="fai", help="fasta index file (.fai) to provide chromosome lengths")
    parser.add_argument("--vcf_path", "-v",dest="vcf_path",help="path to vcf files to use as input, all files ending in .vcf in this path will be parsed")
    parser.add_argument("--out","-o",dest="out_file",help="file to write to for circos R helper script")
    parser.add_argument("--out_detailed","-d",dest="detailed_out_file",help="file to write sample details for events to (for reference only)")
    args = parser.parse_args()
    vcf_path = args.vcf_path
    fai = args.fai

    bin_size=5000000
    bins = createGenomeBins(fai,bin_size)
    bin_break_counts = {} #store a count for every bin
    bin_break_pairs = {} #store the pair details for every event, organized by bin for retrieval
    out_h = open(args.out_file,"w")
    out_d = open(args.detailed_out_file,"w")
    os.chdir(vcf_path)
    sample = ""
    for file in glob.glob("*.vcf"):
        sample = file.rstrip(".vcf")
        print "%s for %s"% (sample,file)
        file = vcf_path + file
        print(file)
        tallyBins(file,bins,bin_break_counts,bin_break_pairs,sample)
    #print bin_break_pairs.keys()
    samp_events = {} #key on sample_id and chr1Bin1chr2Bin2 to collapse any redundant events in the same patient
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
    
                        for linkage in bin_break_pairs[chrom][bin]:
                            (c1,c2,b1,b2,sample,type) = linkage
    
                            if dup_count.has_key(sample):
                                continue
                            dup_count[file] = 1
                            
                            s1 = bins[c1][0][b1]
                            e1 = bins[c1][1][b1]
                            s2 = bins[c2][0][b2]
                            e2 = bins[c2][1][b2]
                            out_string = "chr%s\tchr%s\t%s\t%s\t%s" % (c1,c2,type,s1,s2) 
                            #out_h.write(out_string)
                            if samp_events.has_key(out_string):
                                if samp_events[out_string].has_key(sample):
                                    samp_events[out_string][sample]+=1
                                else:
                                    samp_events[out_string][sample] = 1
                            else:
                                samp_events[out_string] = {sample:1}
            else:
                print "no chromosome for %s" % chrom
    for event in samp_events.keys():
        samps = samp_events[event].keys()
        
        n = len(samps)
        out_string = "%s\t%s\n" % (event,n)
        out_h.write(out_string)
        for samp in samps:
            out_string = "%s\t%s\n" % (samp,event)
            out_d.write(out_string)
def getBin(starts,ends,pos):
    for i in range(0,len(ends)):
        start = starts[i]
        end = ends[i]
        if pos >= start and pos <= end:
            return(i)

    print "%s NOT FOUND!!!" % pos
    print ends[-1]
    exit()

def tallyBins(vcf_file,bins,tally,links,sample):
    #first load chrom, pos details for all and key on ID
    sv_id = {}
    sv_mates = {}
    vcf_reader = vcf.Reader(open(vcf_file,"r"))
    for record in vcf_reader:
        sv_id[record.ID] = record
        if record.INFO.has_key("MATEID"):
            mate_id = record.INFO["MATEID"][0]
            sv_mates[record.ID] = mate_id

                
    for chrom in bins.keys():
        
        
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
                
                
                bin_inds = getBin(bins[chrom][0],bins[chrom][1],record.POS)
                bin_inde = getBin(bins[chrom][0],bins[chrom][1],inv_end)
                
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
                    continue
                mate_bin_ind = getBin(bins[mate_chr][0],bins[mate_chr][1],mate_record.POS)
                twochrom = [mate_chr,chrom]
                twochrom.sort()
                link = twochrom
        
                if twochrom[0] == mate_chr:
                    link.append(mate_bin_ind)
                    link.append(bin_ind)
                else:
                    link.append(bin_ind)
                    link.append(mate_bin_ind)
                link.append(sample)
                link.append("TRA")
        
                if links.has_key(link[0]):
                    if links[link[0]].has_key(link[2]):
                        links[link[0]][link[2]].append(link)
                    else:
                        links[link[0]][link[2]] = [link]
                else:
                    links[link[0]] = {link[2]:[link]}
                
                if tally.has_key(chrom):
                    if tally[chrom].has_key(bin_ind):
                        tally[chrom][bin_ind] +=1
                        
                    else:
                        tally[chrom][bin_ind] = 1
                else:
                    tally[chrom] = {bin_ind:1}

        
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

main()
