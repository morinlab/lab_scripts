import argparse
import math
import re
from collections import deque

"""
Input: a single file containing the segments from all the samples in IGV 
format. Also requires a BED file containing the exon start and end 
positions to create the marker file. NOTE: assumes that chromosome have 
'chr' prefix in the IGV segment file and the BED file.

Output: GISTIC compatible segmentation file with the start position of the
first segment and the end position of the last segment for every chromosome
is adjusted to be the same for all samples. It will also produce a marker
file where the marker positions are derived from the start and end positions
of the segments across samples as well as the start and end positions of the
exons in the BED file.

"""

def main():
    args = parse_args()
    concat_file = args.concatenated_file
    marker_file = args.marker_file
    out_seg = args.out_seg
    bed_file = args.bed_file

    # Parse segs into data structure
    parsed_seg = parse_segs_from_concatenated_file(concat_file)

    # Adjust start and end
    parsed_seg = add_buffer_to_first_seg(parsed_seg)
    parsed_seg = add_buffer_to_last_seg(parsed_seg)

    chrm_min_max = find_chrm_min_max(parsed_seg)

    parsed_seg = adjust_chrm_segs(parsed_seg, chrm_min_max)

    markers = make_markers_from_bed_file(bed_file)
    markers = add_start_and_end_markers(parsed_seg, markers)

    num_marker_dict = find_num_markers(parsed_seg, markers)

    write_marker_file(markers, marker_file)
    write_segmentation_file(parsed_seg, num_marker_dict, out_seg)

    return

def adjust_chrm_segs(parsed_seg, chrm_min_max):
    for sn in parsed_seg.keys():
        for chrm_k in make_chrm_list():
            if chrm_k not in parsed_seg[sn].keys():
                parsed_seg[sn][chrm_k] = deque()
                parsed_seg[sn][chrm_k].append([chrm_min_max[chrm_k][0], chrm_min_max[chrm_k][1],2])
            else:
                parsed_seg[sn][chrm_k][0][0] = chrm_min_max[chrm_k][0]
                parsed_seg[sn][chrm_k][-1][1] = chrm_min_max[chrm_k][1]

    return parsed_seg

def find_chrm_min_max(parsed_seg):
    chrm_min_max = { k : [float("inf"), 0] for k in make_chrm_list() }

    for sn in parsed_seg.keys():
        chrm_dict = parsed_seg[sn]

        for chrm_k in make_chrm_list():

            if chrm_k not in parsed_seg[sn]:
                continue

            chrm_segs = parsed_seg[sn][chrm_k]

            first_seg_start = int(chrm_segs[0][0])
            last_seg_end = int(chrm_segs[-1][1])

            if first_seg_start < chrm_min_max[chrm_k][0]:
                chrm_min_max[chrm_k][0] = first_seg_start

            if last_seg_end > chrm_min_max[chrm_k][1]:
                chrm_min_max[chrm_k][1] = last_seg_end

    for chrm in make_chrm_list():
        if chrm not in chrm_min_max:
            chrm_min_max[chrm] = [1, CHRM_LENGTH_DICT[chrm]]

    return chrm_min_max

def add_buffer_to_first_seg(parsed_seg):
    buffer = 100
    buffer_length = 50

    for sn in parsed_seg:
        chrm_dict = parsed_seg[sn]

        for chrm in make_chrm_list():
            if chrm in chrm_dict:
                seg_start = int(chrm_dict[chrm][0][0])
                chrm_dict[chrm].appendleft([seg_start - buffer, seg_start - buffer + buffer_length, 2])

    return parsed_seg

def add_buffer_to_last_seg(parsed_seg):
    buffer = 100
    buffer_length = 50

    for sn in parsed_seg:
        chrm_dict = parsed_seg[sn]

        for chrm in make_chrm_list():
            if chrm in chrm_dict:
                seg_end = int(chrm_dict[chrm][-1][1])
                chrm_dict[chrm].append([seg_end + buffer_length, seg_end + buffer, 2])

    return parsed_seg

def parse_segs_from_concatenated_file(concat_file):
    parsed_seg = {}

    header_found = False

    p = re.compile("start", re.I)

    for seg in concat_file:

        if not header_found:
            header_found = True
            if re.search(p, seg):
                continue

        fields = seg.rstrip().split('\t')

        sample_name, chrm, start, end = fields[0:4]
        ratio = fields[-1]

        if sample_name not in parsed_seg:
            parsed_seg[sample_name] = {}

        if chrm.startswith('"') and chrm.endswith('"'):
            chrm = chrm[1:-1]

        if chrm not in parsed_seg[sample_name]:
            parsed_seg[sample_name][chrm] = deque()

        parsed_seg[sample_name][chrm].append([start, end, ratio])

    return parsed_seg

def write_segmentation_file(parsed_seg, num_marker_dict, seg_outfile):
    for sn in parsed_seg.keys():
        chrm_dict = parsed_seg[sn]
        for chrm_k in make_chrm_list():
            chrm_segs = parsed_seg[sn][chrm_k]
            num_marker_list = num_marker_dict[sn][chrm_k]
            i = 0
            for chrm_seg in chrm_segs:
                num_marker = num_marker_list[i]
                seg_start, seg_end, logr = chrm_seg[0:3]
                out_string = "{}\t{}\t{}\t{}\t{}\t{}\n".format(sn, chrm_k, seg_start, seg_end, num_marker, logr)
                seg_outfile.write(out_string)
                i += 1
    else:
        seg_outfile.close()

    return

def write_marker_file(markers, marker_outfile):
    counter = 0
    for chrm_k in make_chrm_list():
        chrm_marker_list = markers[chrm_k]
        for pos in chrm_marker_list:
            out_string = "MK_{}\t{}\t{}\n".format(counter, chrm_k, pos)
            counter += 1
            marker_outfile.write(out_string)
    else:
        marker_outfile.close()

    return

def find_num_markers(parsed_seg, markers):
    num_marker_dict = { k : {} for k in parsed_seg.keys() }

    for k in num_marker_dict.keys():
        num_marker_dict[k] = { c : [] for c in make_chrm_list() }

    for sn in parsed_seg.keys():
        chrm_dict = parsed_seg[sn]
        for chrm_k in make_chrm_list():
            chrm_markers = markers[chrm_k]
            segs = chrm_dict[chrm_k]
            for seg in segs:
                marker_count = 0
                for marker in chrm_markers:
                    if marker >= seg[0] and marker <= seg[1]:
                        marker_count += 1
                num_marker_dict[sn][chrm_k].append(marker_count)

    return num_marker_dict

def make_markers_from_bed_file(bed_file):
    bed_markers = { chrm : set() for chrm in make_chrm_list() }

    for line in bed_file:
        chrm, start, end = line.split('\t')[0:3]

        if chrm not in make_chrm_list():
            continue

        bed_markers[chrm].add(start)
        bed_markers[chrm].add(end)

    return bed_markers

def add_start_and_end_markers(parsed_seg, markers):
    out_markers = {}

    for sn in parsed_seg:
        chrm_dict = parsed_seg[sn]
        for chrm_k in make_chrm_list(): 
            chrm_segs = parsed_seg[sn][chrm_k]
            for chrm_seg in chrm_segs:
                markers[chrm_k].add(chrm_seg[0])
                markers[chrm_k].add(chrm_seg[1])

    for chrm_k in make_chrm_list():
        out_markers[chrm_k] = sorted(list(markers[chrm_k]), key=int)

    return out_markers

def make_chrm_list():
    chrm_list = ["chr" + str(i) for i in range(1, 23)]
    return chrm_list

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--concatenated_file', nargs='?', type=argparse.FileType('r'),
                        help='Concatenated IGV formatted sample file.')
    parser.add_argument('-m', '--marker_file', required=True,
                        type=argparse.FileType('w'), help='Specify marker file output.')
    parser.add_argument('-o', '--out_seg', required=True,
                        type=argparse.FileType('w'), help='Specify GISTIC seg file output.')
    parser.add_argument('-b', '--bed_file', required=True,
                        type=argparse.FileType('r'), help='Specify BED file for marker output.')
    args = parser.parse_args()
    return args

CHRM_LENGTH_DICT={
    "chr1":[1, 249250621],
    "chr2":[1, 243199373],
    "chr3":[1, 198022430],
    "chr4":[1, 191154276],
    "chr5":[1, 180915260],
    "chr6":[1, 171115067],
    "chr7":[1, 159138663],
    "chr8":[1, 146364022],
    "chr9":[1, 141213431],
    "chr10":[1, 135534747],
    "chr11":[1, 135006516],
    "chr12":[1, 133851895],
    "chr13":[1, 115169878],
    "chr14":[1, 107349540],
    "chr15":[1, 102531392],
    "chr16":[1, 90354753],
    "chr17":[1, 81195210],
    "chr18":[1, 78077248],
    "chr19":[1, 59128983],
    "chr20":[1, 63025520],
    "chr21":[1, 48129895],
    "chr22":[1, 51304566]
}

if __name__ == '__main__':
    main()
