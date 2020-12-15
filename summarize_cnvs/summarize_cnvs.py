#!/usr/bin/env python


import argparse
import os
import sys
import bisect
import logging
from collections import Counter
from sortedcontainers import SortedSet


class Gene:
    """
    A simple class storing the chromosome, start, end, and name of a gene
    """

    __slots__ = ["name", "chrom", "start", "end"]

    def __init__(self, chrom, start, end, name):
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end


class Chromosome:
    """
    A simple class storying the genomic range of each chromosome, as well as the coordiates of the p and q arm
    """

    def __init__(self, chrom):
        self.chrom = chrom
        self.p_start = None
        self.p_end = None
        self.q_start = None
        self.q_end = None
        self.q_length = 0.0
        self.p_length = 0.0

    def add(self, start, end, arm):

        if arm == "p":
            # Check to see if we already have genomic coordinates for this arm
            if self.p_start is not None or self.p_end is not None:
                raise AttributeError("Coordinates already assigned for %s arm of chromosome \'%s\'" % (arm, self.chrom))
            self.p_start = start
            self.p_end = end
            self.p_length = end - start

            # Sanity check this does not overlap the q arm
            if self.q_start is not None and self.p_end > self.q_start:
                raise AttributeError("For chromosome \'%s\'. The q arm starts before the end of the p arm" % self.chrom)

        elif arm == "q":
            # Check to see if we already have genomic coordinates for this arm
            if self.q_start is not None or self.q_end is not None:
                raise AttributeError("Coordinates already assigned for %s arm of chromosome \'%s\'" % (arm, self.chrom))
            self.q_start = start
            self.q_end = end
            self.q_length = end - start

            # Sanity check this does not overlap the p arm
            if self.p_end is not None and self.q_start < self.p_end:
                raise AttributeError("For chromosome \'%s\'. The q arm starts before the end of the p arm" % self.chrom)
        else:
            # Not the p or q arm. Invalid
            raise AttributeError("Unknown arm specified for chromosome \'%s\':\'%s\'" % (self.chrom, arm))


class SampleCNVs:
    """
    Stores all the copy number events within a given class
    I was going to write this in a function, but having a class is a lot cleaner
    """

    def __init__(self):
        self.starts = {}
        self.ends = {}
        self.cn_states = {}
        self.len_seg = {}
        self.is_chr_prefixed = None
        self.ploidy = None

    def merge_overlapping_cnvs(self, existStarts: iter, existEnds: iter, existCNs: iter, newStart: int, newEnd: int, newCN:int):
        """
        Handle overlapping copy number segments, merging and editing existing segments as necessary

        In essence, if this segment overlaps an existing segment with the same copy number state, those segments are merged
        together. If the copy number state is different, the largest increase is kept,
         truncating or even breaking apart the copy-neutral segment in the process. In the worse case senario,
        this will lead to three different segments being created

        :param oldStart: The start position of the existing segment
        :param oldEnd: The end position of the existing segment
        :param oldCN: The copy number state of the existing segment
        :param newStart: The start position of the new segment
        :param newEnd: The end position of the new segment
        :param newCN: The copy number state of the new segment
        :return: A tuple containing the new and old events. Ex {[newEvent1, newEvent2], [oldEvent1, oldEvent2]}
        """

        def reduce_new_event(existStarts, existEnds, existCNs, newStart, newEnd, newCN, oldStart, oldEnd):

            # The existing event has higher priority than the old event
            # Truncate the newer event
            # Since this could break the newer event apart into two pieces, lets do that now, and see
            # if the output segments are valid
            outStart1 = newStart
            outEnd1 = oldStart
            outStart2 = oldEnd
            outEnd2 = newEnd
            if outStart1 < outEnd1:  # Valid segment, process
                existStarts, existEnds, existCNs = self.merge_overlapping_cnvs(existStarts, existEnds, existCNs,
                                                                             outStart1, outEnd1, newCN)
            if outStart2 < outEnd2:  # Valid segment
                existStarts, existEnds, existCNs = self.merge_overlapping_cnvs(existStarts, existEnds, existCNs,
                                                                             outStart2, outEnd2, newCN)
            return existStarts, existEnds, existCNs

        def reduce_old_event(existStarts, existEnds, existCNs, newStart, newEnd, newCN, oldStart, oldEnd, oldCN):

            # The newer event is a "higher"-level deletion
            # Split/Truncate the older event
            outStart1 = oldStart
            outEnd1 = newStart
            outStart2 = newEnd
            outEnd2 = oldEnd
            outCN = oldCN
            # Delete existing event
            existStarts.pop(start_bisect)
            existEnds.pop(start_bisect)
            existCNs.pop(start_bisect)

            if outStart2 < outEnd2:  # Valid segment, do the last one first to maintain sorted order
                existStarts.insert(start_bisect, outStart2)
                existEnds.insert(start_bisect, outEnd2)
                existCNs.insert(start_bisect, outCN)
            if outStart1 < outEnd1:  # Also valid
                existStarts.insert(start_bisect, outStart1)
                existEnds.insert(start_bisect, outEnd1)
                existCNs.insert(start_bisect, outCN)

            # Check for any more overlaps
            return self.merge_overlapping_cnvs(existStarts, existEnds, existCNs, newStart, newEnd, newCN)

        # First, does this new segment actually overlap anything?
        start_bisect = bisect.bisect_right(existEnds, newStart)
        end_bisect = bisect.bisect_left(existStarts, newEnd)

        if start_bisect == end_bisect:
            # There are no overlaps. Simply add this new event in, and return it
            existStarts.insert(start_bisect, newStart)
            existEnds.insert(start_bisect, newEnd)
            existCNs.insert(start_bisect, newCN)
            return existStarts, existEnds, existCNs

        # Grab the first overlap
        oldStart = existStarts[start_bisect]
        oldEnd = existEnds[start_bisect]
        oldCN = existCNs[start_bisect]

        # Simplest case first. If both these events have the same CN state, lets merge them together
        if oldCN == newCN:

            outStart = oldStart if oldStart < newStart else newStart
            outEnd = oldEnd if oldEnd > newEnd else newEnd
            # Delete the existing event
            existStarts.pop(start_bisect)
            existEnds.pop(start_bisect)
            existCNs.pop(start_bisect)

            # Check for any more overlaps
            return self.merge_overlapping_cnvs(existStarts, existEnds, existCNs, outStart, outEnd, newCN)
        else:
            # These segments overlap and have different CN states.
            # Lets keep the highest-level event
            if newCN <= 2 and oldCN <= 2: # Deletion
                if oldCN < newCN:
                    # The older event is a "higher"-level deletion
                    return reduce_new_event(existStarts, existEnds, existCNs, newStart, newEnd, newCN, oldStart, oldEnd)
                else:
                    return reduce_old_event(existStarts, existEnds, existCNs, newStart, newEnd, newCN, oldStart, oldEnd, oldCN)
            if newCN >= 2 and oldCN >= 2:  # Gain/Amp
                if oldCN > newCN:
                    # The older event is a higher level gain. Split/truncate the new event
                    return reduce_new_event(existStarts, existEnds, existCNs, newStart, newEnd, newCN, oldStart,
                                            oldEnd)
                else:
                    return reduce_old_event(existStarts, existEnds, existCNs, newStart, newEnd, newCN, oldStart, oldEnd,
                                          oldCN)
            else:
                # One event must be a gain/amp, while the other is a deletion. In this case, keep both, and subset the
                # larger event by the smaller one
                if oldEnd - oldStart < newEnd - newStart:
                    # The older event is smaller. Split the newer event
                    return reduce_new_event(existStarts, existEnds, existCNs, newStart, newEnd, newCN, oldStart,
                                            oldEnd)
                else:
                    # The newer event is smaller. We should split the older event
                    return reduce_old_event(existStarts, existEnds, existCNs, newStart, newEnd, newCN, oldStart, oldEnd,
                                          oldCN)


    def add(self, chrom: str, start: int, end: int, cn: float):

        if isinstance(cn, float):
            cn = round(cn)

        assert start < end
        # Have we seen any events for this chromosome before?
        if chrom not in self.len_seg:
            # If not, just store this segment. Simple!
            self.starts[chrom] = [start]
            self.ends[chrom] = [end]
            self.cn_states[chrom] = [cn]
            self.len_seg[chrom] = 1

            # Are we using chr-prefixed contig names?
            if self.is_chr_prefixed is None:
                self.is_chr_prefixed = True if chrom.startswith("chr") else False
        else:
            # If we have seen events on this chromosome before, we need to compare this new segment with existing segments
            # In theory, events should be non-overlapping, but I don't want to assume that because then things will break
            # horribly later

            self.starts[chrom], self.ends[chrom], self.cn_states[chrom] = \
                self.merge_overlapping_cnvs(self.starts[chrom], self.ends[chrom], self.cn_states[chrom], start, end, cn)
            self.len_seg[chrom] = len(self.cn_states)


    def overlap_chrom(self, chromosome: Chromosome, threshold: float = 0.8):

        """
        Calculates the overlap between the segments stored here and a given chromosome
        :param chromosome: A Chromosome() object defining the start and end of the p and q arms, respectively
        :return: A dictionary containing
        """

        # Do the CNVs within this chromosome encompass a given chromosome more than the specified threshold?
        chrom_name = chromosome.chrom

        # if the p or q arm coordinates are not set, use a placeholder
        if chromosome.p_start is None:
            chromosome.p_start = -100000
            chromosome.p_end = -100000
            chromosome.p_length = 1
        if chromosome.q_start is None:
            chromosome.q_start = -100000
            chromosome.q_end = -100000
            chromosome.q_length = 1

        # Handle "chr" prefix nonsense
        if self.is_chr_prefixed and not chrom_name.startswith("chr"):
            chrom_name = "chr" + chrom_name
        elif not self.is_chr_prefixed and chrom_name.startswith("chr"):
            chrom_name = chrom_name.replace("chr", "")

        try:
            chrom_starts = self.starts[chrom_name]
        except KeyError:
            # No CN events were provided for this chromosome, hence there can be no overlap
            return {}
        chrom_ends = self.ends[chrom_name]
        chrom_cn = self.cn_states[chrom_name]

        homdel_sum = {"p": 0, "q": 0, "chrom": 0}
        del_sum = {"p": 0, "q": 0, "chrom": 0}
        gain_sum = {"p": 0, "q": 0, "chrom": 0}
        amp_sum = {"p": 0, "q": 0, "chrom": 0}

        # Check for overlapping segments for each arm
        for (start, end, cn) in zip(chrom_starts, chrom_ends, chrom_cn):

            # Ignore copy-neutral segments
            if cn == 2:
                continue

            if end < chromosome.p_start or start > chromosome.q_end:
                # Segment falls out of range
                continue
            if start < chromosome.p_end:
                olap_start = start if start > chromosome.p_start else chromosome.p_start
                olap_end = end if end < chromosome.p_end else chromosome.p_end
                if cn < 2:
                    del_sum["p"] += olap_end - olap_start
                    del_sum["chrom"] +=  olap_end - olap_start
                    if cn < 1:
                        homdel_sum["p"] += olap_end - olap_start
                        homdel_sum["chrom"] += olap_end - olap_start
                elif cn > 2:
                    gain_sum["p"] += olap_end - olap_start
                    gain_sum["chrom"] += olap_end - olap_start
                    if cn > 3:
                        amp_sum["p"] += olap_end - olap_start
                        # amp_sum["chrom"] += olap_end - olap_start

            if end > chromosome.q_start:  # We use an if, not elif, in case a segment overlaps both the p and q arm
                olap_start = start if start > chromosome.q_start else chromosome.q_start
                olap_end = end if end < chromosome.q_end else chromosome.q_end
                if cn < 2:
                    del_sum["q"] += olap_end - olap_start
                    del_sum["chrom"] +=  olap_end - olap_start
                    if cn < 1:
                        homdel_sum["q"] += olap_end - olap_start
                        homdel_sum["chrom"] += olap_end - olap_start
                elif cn > 2:
                    gain_sum["q"] += olap_end - olap_start
                    gain_sum["chrom"] += olap_end - olap_start
                    if cn > 3:
                        amp_sum["q"] += olap_end - olap_start
                        # amp_sum["chrom"] += olap_end - olap_start

        events = {}
        # Now, calculate the fraction of each overlap
        # Start from the highest level and biggest event, and work our way down/smaller
        # Amplifications
        wholeChromAmp = False
        if amp_sum["p"] / chromosome.p_length > threshold:  # p Amp
            if amp_sum["q"] / chromosome.q_length > threshold:  # Whole chromosome Amp
                events[chrom_name + "chrom"] = "AMP"
                wholeChromAmp = True
            events[chrom_name + "p"] = "AMP"
        if amp_sum["q"] / chromosome.q_length > threshold:  # q Amp
            events[chrom_name + "q"] = "AMP"
        # Gains
        if not wholeChromAmp:
            if chrom_name + "p" not in events and gain_sum["p"] / chromosome.p_length > threshold:  # p Gain
                if gain_sum["q"] / chromosome.q_length > threshold:  # Whole chromosome gains
                    events[chrom_name + "chrom"] = "GAIN"
                events[chrom_name + "p"] = "GAIN"
            if chrom_name + "q" not in events and gain_sum["q"] / chromosome.q_length > threshold:  # q Gain
                events[chrom_name + "q"] = "GAIN"

        # Homozygous deletions
        wholeChromDel = False
        if homdel_sum["p"] / chromosome.p_length > threshold:  # p HOMDEL
            if homdel_sum["q"] / chromosome.q_length > threshold:  # Whole chromosome homdel
                events[chrom_name + "chrom"] = "HOMDEL"
                wholeChromDel = True
            events[chrom_name + "p"] = "HOMDEL"
        if homdel_sum["q"] / chromosome.q_length > threshold:  # q HOMDEL
            events[chrom_name + "q"] = "HOMDEL"
        # Heterozygous deletions
        if not wholeChromDel:
            if chrom_name + "p" not in events and del_sum["p"] / chromosome.p_length > threshold:  # p deletion
                if del_sum["q"] / chromosome.q_length > threshold:  # Whole chromosome hetdel
                    events[chrom_name + "chrom"] = "HETLOSS"
                events[chrom_name + "p"] = "HETLOSS"
            if chrom_name + "q" not in events and del_sum["q"] / chromosome.q_length > threshold:  # q deletions
                events[chrom_name + "q"] = "HETLOSS"

        return events

    def get_overlap_regions(self, gene_coords: dict, size_threshold = 20000000, olap_threshold = 0.6):

        olap_genes = {}
        for chrom in self.starts.keys():
            # Handle chr-prefix shenanigans
            g_chrom = chrom
            is_chr_prefixed = next(iter(gene_coords.keys())).startswith("chr")
            if is_chr_prefixed and not self.is_chr_prefixed:
                g_chrom = "chr" + chrom
            elif not is_chr_prefixed and self.is_chr_prefixed:
                g_chrom = chrom.replace("chr", "")

            try:
                genes_on_chrom = gene_coords[g_chrom]
            except KeyError:
                # No genes are on this contig
                continue

            olap_genes_chrom = {}
            for start, end, cn_state in zip(self.starts[chrom], self.ends[chrom], self.cn_states[chrom]):

                # Ignore copy-neutral events
                if cn_state == 2:
                    continue
                # Ignore events which are larger than the specified threshold
                if end - start > size_threshold:
                    continue

                for gene_name, coords in genes_on_chrom.items():
                    # Does this gene/region overlap this CNV?
                    if start > coords.end or end < coords.start:
                        continue  # Not overlapping
                    # Calculate the amount of overlap
                    olap_start = start if start > coords.start else coords.start
                    olap_end = end if end < coords.end else coords.end
                    # Store this overlap
                    # Here we are assuming no overlapping segments, since those were delt with earlier
                    if gene_name not in olap_genes_chrom:
                        olap_genes_chrom[gene_name] = {}
                    if cn_state not in olap_genes_chrom[gene_name]:
                        olap_genes_chrom[gene_name][cn_state] = olap_end - olap_start
                    else:
                        olap_genes_chrom[gene_name][cn_state] += olap_end - olap_start

            # Now that we have processed all CNVs, see which events are significantly overlapping
            for gene_name, olap_by_cn in olap_genes_chrom.items():
                # Get the size of this gene
                gene_size = gene_coords[g_chrom][gene_name].end - gene_coords[g_chrom][gene_name].start

                # Is this gene gained?
                amp_size = sum(olap_by_cn[x] for x in olap_by_cn.keys() if x > 3)
                if amp_size / gene_size > olap_threshold:
                    # Gene is amplified
                    olap_genes[gene_name] = "AMP"
                else:
                    try:
                        gain_size = amp_size + olap_by_cn[3]
                        if gain_size / gene_size > olap_threshold:
                            # Gene is gained
                            olap_genes[gene_name] = "GAIN"
                    except KeyError:  # i.e. this gene is not gained
                        pass

                # Is this gene deleted?
                # NOTE: THIS WILL BREAK IF olap_threshold IS LESS THAN 0.5, AS THEN A GENE CAN BE BOTH GAINED AND DELETED
                # I am using a dictionary, so the deletions will override gains
                homdel_size = sum(olap_by_cn[x] for x in olap_by_cn.keys() if x < 1)
                if homdel_size / gene_size > olap_threshold:
                    # Gene has a homozygous deletion
                    olap_genes[gene_name] = "HOMDEL"
                else:
                    try:
                        hetdel_size = homdel_size + olap_by_cn[1]
                        if hetdel_size / gene_size > olap_threshold:
                            # Gene has a heterozygous deletion
                            olap_genes[gene_name] = "HETLOSS"
                    except KeyError:  # i.e. this gene does not have a heterozygous deletion
                        pass

        return olap_genes


    def adjust_ploidy(self, arm_coords, sample_name, redo=False):
        """
        Adjust a sample's CN states based upon ploidy

        Calculate the average copy number state across the entire genome. If a genome has a higher ploidy state, all the
        copy number values are adjusted based upon the average.
        :param arm_coords: A dictionary containing {"chromosome_name": Chromosome()} objects which list chromosomal coordinates
        :param sample_name: A string specifying the sample name. Used for debugging/error purposes
        :param redo: Should we override any existing ploidy state for this sample?
        :return: None
        """

        if self.ploidy is not None and not redo:
            # Ploidy for this sample has already been calculated. Don't do anything
            return None

        # Store the length affected for each ploidy state
        ploidy_cov = {}
        genome_size = 0

        for chromosome in arm_coords.values():
            # Find overlapping CNVs for this arm

            # if the p or q arm coordinates are not set, use a placeholder
            if chromosome.p_start is None:
                chromosome.p_start = -100000
                chromosome.p_end = -100000
                chromosome.p_length = 1
            if chromosome.q_start is None:
                chromosome.q_start = -100000
                chromosome.q_end = -100000
                chromosome.q_length = 1

            # Save the size of this chromosome
            genome_size += chromosome.p_length
            genome_size += chromosome.q_length

            chrom_name = chromosome.chrom
            # Handle "chr" prefix nonsense
            if self.is_chr_prefixed and not chrom_name.startswith("chr"):
                chrom_name = "chr" + chrom_name
            elif not self.is_chr_prefixed and chrom_name.startswith("chr"):
                chrom_name = chrom_name.replace("chr", "")
            try:
                chrom_starts = self.starts[chrom_name]
            except KeyError:
                # No CN events were provided for this chromosome, hence there can be no overlap
                continue
            chrom_ends = self.ends[chrom_name]
            chrom_cn = self.cn_states[chrom_name]

            # Check for overlapping segments for each arm
            for (start, end, cn) in zip(chrom_starts, chrom_ends, chrom_cn):

                # Check p-arm overlap
                if end < chromosome.p_start or start > chromosome.q_end:
                    # Segment falls out of range
                    continue
                if start < chromosome.p_end:
                    olap_start = start if start > chromosome.p_start else chromosome.p_start
                    olap_end = end if end < chromosome.p_end else chromosome.p_end

                    # Store the length of this overlap and associated CN
                    if cn not in ploidy_cov:
                        ploidy_cov[cn] = 0
                    ploidy_cov[cn] += olap_end - olap_start

                if end > chromosome.q_start:  # We use an if, not elif, in case a segment overlaps both the p and q arm
                    olap_start = start if start > chromosome.q_start else chromosome.q_start

                    olap_end = end if end < chromosome.q_end else chromosome.q_end
                    # Store the length of this overlap and associated CN
                    if cn not in ploidy_cov:
                        ploidy_cov[cn] = 0
                    ploidy_cov[cn] += olap_end - olap_start

        # Now that we have calculated the number of bases affected by each CN state, calculate the ploidy
        x = 0
        for cn, bases in ploidy_cov.items():
            x += cn * bases

        ploidy = x / genome_size
        av_ploidy = round(ploidy)

        # Sanity check
        if av_ploidy < 1:
            raise TypeError("Ploidy of sample \'%s\' was calculated to be below 1!" % sample_name)

        # If this tumour is not diploid, adjust the ploidy
        if av_ploidy != 2:
#            logging.info("\'%s\' was calculated to have a ploidy of %s" % (sample_name, av_ploidy))
            print("\'%s\' was calculated to have a ploidy of %s" % (sample_name, av_ploidy))
            ploidy_dif = av_ploidy - 2
            new_cns = {}
            for chrom, cn in self.cn_states.items():
                x = list(y - ploidy_dif for y in cn)
                new_cns[chrom] = x
            self.cn_states = new_cns

        self.ploidy = av_ploidy


def get_args():

    def is_valid_file(path, parser):
        """
        Checks to ensure the specified file exists
        :param path: A string containing a filepath
        :param parser: An argparse.ArgumentParser() object
        :return: path, if the string is a valid directory
        :raises: parser.error() if the file does not exist
        """
        if os.path.exists(path) and os.path.isfile(path):
            return path
        else:
            raise parser.error("Unable to locate \'%s\': No such file or directory" % path)

    epilog = os.linesep.join(["The --cnvs file should contain the following colummns: Tumor_Sample_Barcode, chromosome, start, end, CN.",
                              "The --arms file should contain the following columns: chromosome, start, end, arm.",
                              "If a gene is partially overlapping a deletion, that gene is considered to be deleted, while genes partially overlapping gains are not considered gained."
                              ])
    parser = argparse.ArgumentParser(description="Summarizes CNVs on a per-gene and chromosomal basis. Handles overlapping events, and adjusts for ploidy", epilog=epilog)

    input = parser.add_argument_group("Input files")
    input.add_argument("-c", "--cnvs", metavar="SEG", default=None, required=True, type=lambda x: is_valid_file(x, parser), help="Input tab-delimited file summarizing copy number events (with absolute CN status)")
    input.add_argument("-g", "--genes", metavar="BED", default=None, required=True, type=lambda x: is_valid_file(x, parser), help="Input BED4+ file listing start and end positions of genes/exons")
    input.add_argument("-a", "--arms", metavar="TSV", default=None, required=True, type=lambda x: is_valid_file(x, parser), help="Input tab-delimited file listing the positions of chromosome arms")

    parser.add_argument("-o", "--output", metavar="TSV", default=None, required=True, type=str, help="Output file listing all genes and chromosomal arms affected by CNVs")
    parser.add_argument("--olap_threshold", metavar=0.7, default=0.7, type=float, help="--genes and --arms must have this fraction of their length affected by CNVs before they are called as deleted/gained")
    parser.add_argument("--focal_threshold", metavar="INT", default=30000000, type=int, help="CNVs overlapping --genes will only be counted if they are smaller than this [Default:30000000]")

    args = parser.parse_args()

    if args.olap_threshold > 1.0:
        raise parser.error("\'--olap_threshold\' can not be greater than 1")
    return args

def load_gene_coords_bed(bed_file):
    """
    Load in the genomic coordinates of genes in the human genome from the specified BED4+ file.

    The following columns are required (no header)
    chrom   start   end gene

    A gene can be split across multiple BED entries (ex. exonic coordinates). In this case, the start and end of the first and last exon,
    respectively, will be used as the gene start/end coordinates

    :param bed_file: A string containing a filepath to a BED4+ file listing the positions of genes
    :return: A dictionary storing {gene_name: Gene()}
    """

    gene_coords = {}
    skipped_genes = []

    i = 0
    with open(bed_file) as f:
        for line in f:
            i += 1
            line = line.rstrip("\n").rstrip("\r")  # Remove line endings
            cols = line.split("\t")

            # Skip empty lines
            if not line:
                continue

            # Since this BED file could contain more than 4 columns (ex. a BED6 or BED12 file), manually unpack and inspect the first
            # four columns, since those are the columns we care about
            try:
                chrom = cols[0]
                start = cols[1]
                end = cols[2]
                gene = cols[3]
            except IndexError as e:
                # This BED entry was trucated
                raise AttributeError("Unable to parse line %s of \'%s\' as it contains less than four columns (chrom, start, end, gene)" % (i, bed_file)) from e

            # Check that the start and end are actually genomic coordinates
            if not start.isdigit():
                raise TypeError("When parsing line %s of %s, the start position \'%s\' is not a valid genomic coordinate" % (i, bed_file, start))
            if not end.isdigit():
                raise TypeError("When parsing line %s of %s, the end position \'%s\' is not a valid genomic coordinate" % (i, bed_file, end))

            start = int(start)
            end = int(end)

            # Have we seen this chromosome before?
            if chrom not in gene_coords:
                gene_coords[chrom] = {}
            # Have we seen this gene before?
            if gene in gene_coords[chrom]:
                # If so, then we need to update the start/end of this gene based on this new BED entry
                existing_gene_entry = gene_coords[chrom][gene]

                # Sanity check
                if chrom != existing_gene_entry.chrom:
                    raise AttributeError("We found two entries for gene \'%s\'. One is found on chromosome \'%s\', while the other is found on \'%s\'"
                                         % (gene, chrom, existing_gene_entry.chrom))
                if start < existing_gene_entry.start:
                    existing_gene_entry.start = start
                if end > existing_gene_entry.end:
                    existing_gene_entry.end = end
            else:
                # If we haven't seen this gene before, create a new entry for it
                gene_coords[chrom][gene] = Gene(chrom, start, end, gene)

    # Now that we have processed all genes, make sure we actually found Entrez IDs for a handful of genes
    # If we haven't, it means the input file likely didn't contain Hugo Symbols
    if len(gene_coords) == 0:
        raise AttributeError("No Hugo Symbols from the --genes file (column 4) were found in the --entrez_ids file. An example gene is %s" % (skipped_genes[0]))
    elif len(skipped_genes) > 0:
        skipped_genes = set(skipped_genes)  # Remove duplicate entries (ex. multiple exons corresponding to the same gene)
        sys.stderr.write("WARNING: %s genes in the --genes file did not have a corresonding Entrez ID" % len(skipped_genes))
        sys.stderr.write(os.linesep)

    return gene_coords


def load_chrom_arm(arm_file):
    """
    Loads in the genomic coordinates corresponding to the arm of each chromosome

    The input should be a tab-delimited file with the following columns
    chromosome  start   end arm

    :param arm_file: A string containing a filepath to a tab-delimited file containing genomic coordinates for chromosome arms
    :return: A dictionary storing the genomic range of each chromosome, and each individual arm
    """

    arm_chrom = {}
    required_cols = ["chromosome", "start", "end", "arm"]
    header_cols = {}

    i = 0
    with open(arm_file) as f:
        for line in f:
            i += 1
            line = line.rstrip("\n").rstrip("\r")  # Remove line endings
            cols = line.split("\t")

            # Skip empty lines
            if not line:
                continue

            # If we haven't parsed the header yet, assume this is the first line of the file (aka the header)
            if not header_cols:
                j = 0
                for col in cols:
                    if col in required_cols:
                        header_cols[col] = j
                    j += 1

                # Check to make sure all required columns are found
                for col in required_cols:
                    if col not in header_cols:
                        raise AttributeError("Unable to locate column %s in the chromosome arm positions file \'%s\'" % (col, arm_file))
                # If we get this far, the header is valid
                continue

            try:
                arm_attributes = {x: cols[y] for x, y in header_cols.items()}
            except IndexError:
                # We were unable to find a required column. This line is likely truncated
                raise AttributeError("Unable to parse line %s of the chromosome arm positions file %s: The line appears truncated" % (i, arm_file))

            # Have we seen this chromosome before? If not, lets create an object to store its genomic features
            chrom = arm_attributes["chromosome"]
            if chrom not in arm_chrom:
                arm_chrom[chrom] = Chromosome(chrom)

            # Save the arm coordinates in the chromosome
            try:
                arm_chrom[chrom].add(int(arm_attributes["start"]), int(arm_attributes["end"]), arm_attributes["arm"])
            except ValueError as e:
                raise TypeError("Unable to process line %s of \'%s\': start and end must be integers") from e

    return arm_chrom


"""
def get_overlap_genes(chrom: str, start: int, end: int, copy_num: int, gene_cords: dict):

    #""
    Find the genes which overlap a given segment

    This is a very brute-force approach which will check all genes on the target chromosome for overlap. Pre-sorting and bisecting genomic
    regions would be significantly faster, but 1) You need to handle overlapping genes, and 2) I am assuming the performance penalty won't
    matter too much.

    Note that gains and losses are handled differently if a segment partially overlaps a gene. If the event is a loss, the gene is included.
    If the event is a gain, the gene is NOT included, as that copy is likely not functional

    :param chrom: A string coresponding to the contig name
    :param start: An int specifying the start of the segment
    :param end: An int specifying the end of the segment
    :param copy_num: An integer specifying the copy number of this segment
    :param gene_cords: A dictionary storing the positions of genes, in the format {chromosome: {gene1: attr, gene2: attr...}}
    :return:
    #""

    # Handle chr-prefix shenanigans
    is_chr_prefixed = next(iter(gene_cords.keys())).startswith("chr")
    if is_chr_prefixed and not chrom.startswith("chr"):
        chrom = "chr" + chrom
    elif not is_chr_prefixed and chrom.startswith("chr"):
        chrom = chrom.replace("chr", "")

    try:
        genes_on_chrom = gene_cords[chrom]
    except KeyError:
        # No genes are on this contig. This could be a user error, or this a contig with no annotated genes
        return []

    # Go through each gene on this chromosome and see if the coordinates overlap with our regions
    olap_genes = []
    for entrez_id, gene_info in genes_on_chrom.items():
        if start < gene_info.end and end > gene_info.start:
            # Overlap found
            # Now for the tricky part
            # If the segment partially overlaps this gene, then we need to handle gains and losses differently
            # Deleting or gaining half of a gene will likely cause it to no longer function
            # Thus, partial deletions of genes could be drivers, while partial gains of genes are likely never drivers
            if copy_num < 2:
                olap_genes.append(entrez_id)
            elif copy_num > 2:
                if start > gene_info.start or end < gene_info.end:
                    continue  # Partially overlapping
                else:
                    olap_genes.append(entrez_id)
            else:
                raise NotImplementedError("Not calculating overlapping genes for copy-neutral segments")

    return olap_genes
"""

def generate_cnv_files(cnv_segs, gene_regions_bed, arm_regions, outfile, olap_threshold:float = 0.6, focal_cn_thresh:int = 30000000):
    """
    Characterize focal and arm-level copy number events, and summarize them in the respective output files.

    For focal events (i.e. events smaller than the specified threshold), identify all genes which overlap the even, and save them to the out_cnv_gene file
    All events are used to flag chromosomes or individual chromosomal arms which are gained or amplified.

    The cnv_segs file should have the following columns (extra columns are ignored):
    Tumor_Sample_Barcode    chromosome  start"  end CN

    Where CN is the absolute copy number

    The gene_regions_bed can have multiple entries for each gene (ex. exons). The maximum and minimum of the cumulative entries is used to define the
    gene boundaries

    The arm_regions file should contain the following columns:
    chromosome  start   end arm

    :param cnv_segs: A string containing a filepath to a tab-delimited file containing copy number events
    :param gene_regions_bed: A string specifying a filepath to a BED4+ file containing gene positions
    :param arm_regions: A string specifying a filepath to a tab-delimited file specifying the positions of chromosomal arms
    :param outfile: A string specifying the output file name. Each row will list one gene affected by a CNV, as well as chromosomes and chromosomal arms
    :param olap_threshold: A float specifying the proportion of a chromosome/arm/region which has to be affected by CNVs before the region is flagged
    :param focal_cn_thresh: The maximum size of an event for it to be considered "focal", in bp
    :return: A list of samples which have CNV information
    """

    # First things first, lets load in the gene regions file and figure out where each gene is
    gene_coords = load_gene_coords_bed(gene_regions_bed)

    # Load in the chromosome arm regions
    arm_coods = load_chrom_arm(arm_regions)

    # Process copy number segments
    required_cols = ["Tumor_Sample_Barcode", "chromosome", "start", "end", "CN"]
    header_cols = {}

    sample_cnvs = {}

    i = 0
    with open(cnv_segs) as f:
        for line in f:
            i += 1
            line = line.rstrip("\n").rstrip("\r")  # Remove line endings
            cols = line.split("\t")

            # Skip empty lines
            if not line:
                continue

            # If we haven't parsed the header yet, assume this is the first line of the file (aka the header)
            if not header_cols:
                j = 0
                for col in cols:
                    if col in required_cols:
                        header_cols[col] = j
                    j += 1

                # Check to make sure all required columns are found
                for col in required_cols:
                    if col not in header_cols:
                        raise AttributeError("Unable to locate column \'%s\' in the CNV segments file \'%s\'" % (col, cnv_segs))
                # If we get this far, the header is valid
                continue

            # Process this CNV entry
            cnv_attributes = {x: cols[y] for x, y in header_cols.items()}

            # Have we processed CNVs from this sample before?
            if cnv_attributes["Tumor_Sample_Barcode"] not in sample_cnvs:
                sample_cnvs[cnv_attributes["Tumor_Sample_Barcode"]] = SampleCNVs()
            # Sterilize input and ensure it is valid, and store these events
            try:
                cnv_attributes["CN"] = float(cnv_attributes["CN"])
                sample_cnvs[cnv_attributes["Tumor_Sample_Barcode"]].add(
                    cnv_attributes["chromosome"],
                    int(cnv_attributes["start"]),
                    int(cnv_attributes["end"]),
                    cnv_attributes["CN"])
            except ValueError as e:
                raise TypeError("Unable to process line %s of \'%s\': start, end, and CN must be integers or floats" % (i, cnv_segs)) from e


    # Now adjust for ploidy of each case
    for sample, cnvs in sample_cnvs.items():
        cnvs.adjust_ploidy(arm_coods, sample)

    # Now that we have processed all CNVs, lets see which genes have events, and write those out
    with open(outfile, "w") as o:

        # Write output file header
        out_header = ["Sample", "Region", "Event"]
        o.write("\t".join(out_header))
        o.write(os.linesep)

        # Process segments
        for sample, cnvs in sample_cnvs.items():
            olap_genes = cnvs.get_overlap_regions(gene_coords, focal_cn_thresh, olap_threshold)

            for gene_name, cn in olap_genes.items():
                out_line = [sample, gene_name, cn]
                o.write("\t".join(out_line))
                o.write(os.linesep)


        # Now that we have processed all the CNVs, identify which samples have arm-level and whole chromosomal copy number changes
        for sample, cnvs in sample_cnvs.items():
            for chrom, chrom_info in arm_coods.items():

                # For now, skip sex chromosomes
                if chrom == "X" or chrom == "chrX" or chrom == "Y" or chrom == "chrY":
                    continue

                events = cnvs.overlap_chrom(chrom_info, threshold=olap_threshold)
                # If there are large-scale CNVs in this sample, output them to the Arm flat file
                for arm, type in events.items():
                    out_line = [
                        sample,
                        arm,
                        type
                    ]
                    o.write("\t".join(out_line))
                    o.write(os.linesep)

    return list(sample_cnvs.keys())


def main(args=None):

    # Obtain input
    if args is None:
        args = get_args()

    generate_cnv_files(args.cnvs, args.genes, args.arms, args.output, args.olap_threshold, args.focal_threshold)

if __name__ == "__main__":
    main()
