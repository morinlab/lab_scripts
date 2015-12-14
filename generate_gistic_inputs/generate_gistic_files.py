import argparse
import math

"""

Input: a two column, tab delimited file where the first column contains the
absolute paths of the Sequenza seg file and the second column contains the
sample name for the associated Sequenza seg file.
Output: GISTIC compatible segmentation file with the start position of the
first segment and the end position of the last segment for every chromosome
is adjusted to be the same for all samples. It will also produce a marker
file where the markers are the set of start positions of the segments
from all samples. 

TODO: Add functionality that uses the BAF position from Sequenza segmentation 
output as markers instead of segment start and end position.

"""

def main():
    args = parse_args()
    input_table = args.input_table
    marker_file = args.marker_file
    out_seg = args.out_seg

    seg_files = []
    sample_names = []

    for line in input_table:
        file_name, name = line.split('\t')
        file = open(file_name, 'r')
        seg_files.append(file)
        sample_names.append(name.rstrip())

    # TODO: Make sure that samples names are unique.
    if not len(seg_files) == len(sample_names):
        print "The number of sample names do not match the number of specified segmentation files."
        return

    parsed_seg = parse_segs(seg_files, sample_names)

    markers = make_markers(parsed_seg)

    chrm_min_max = find_chrm_min_max(parsed_seg)

    parsed_seg = adjust_chrm_segs(parsed_seg, chrm_min_max) 

    num_marker_dict = find_num_markers(parsed_seg, markers)

    write_marker_file(markers, marker_file)
    write_segmentation_file(parsed_seg, num_marker_dict, out_seg)

    return

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

def make_markers(parsed_seg):
    markers = { k : set() for k in make_chrm_list() }
    out_markers = {}

    for sn in parsed_seg.keys():
        chrm_dict = parsed_seg[sn]
        for chrm_k in make_chrm_list():
            chrm_segs = parsed_seg[sn][chrm_k]
            for chrm_seg in chrm_segs:
                markers[chrm_k].add(chrm_seg[0])
                markers[chrm_k].add(chrm_seg[1])

    for chrm_k in make_chrm_list():
        sorted_marker_list = sorted(list(markers[chrm_k]), key=int)
        out_markers[chrm_k] = sorted_marker_list

    return out_markers

def adjust_chrm_segs(parsed_seg, chrm_min_max):
    for sn in parsed_seg.keys():
        for chrm_k in make_chrm_list():
            parsed_seg[sn][chrm_k][0][0] = chrm_min_max[chrm_k][0]
            parsed_seg[sn][chrm_k][-1][1] = chrm_min_max[chrm_k][1]

    return parsed_seg


def find_chrm_min_max(parsed_seg):
    chrm_min_max = { k : [float("inf"), 0] for k in make_chrm_list() }

    for sn in parsed_seg.keys():
        chrm_dict = parsed_seg[sn]

        for chrm_k in make_chrm_list():
            chrm_segs = parsed_seg[sn][chrm_k]

            first_seg_start = int(chrm_segs[0][0])
            last_seg_end = int(chrm_segs[-1][1])

            if first_seg_start < chrm_min_max[chrm_k][0]:
                chrm_min_max[chrm_k][0] = first_seg_start

            if last_seg_end > chrm_min_max[chrm_k][1]:
                chrm_min_max[chrm_k][1] = last_seg_end

    return chrm_min_max


def parse_segs(seg_files, sample_names):
    parsed_seg = {}
    sample_index = 0

    for seg_file in seg_files:

        header_found = False
        sample_name = sample_names[sample_index]
        sample_index += 1
        parsed_seg[sample_name] = {}

        for seg in seg_file:

            if not header_found:
                header_found = True
                continue

            fields = seg.split('\t')
            chrm, start, end = fields[0:3]
            log_ratio = math.log(float(fields[6]), 2) - 1

            if chrm.startswith('"') and chrm.endswith('"'):
                chrm = chrm[1:-1]

            if chrm not in parsed_seg[sample_name]:
                parsed_seg[sample_name][chrm] = []

            parsed_seg[sample_name][chrm].append([start, end, log_ratio])
    else:
        seg_file.close()

    return parsed_seg

def make_chrm_list():
    chrm_list = [str(i) for i in range(1, 23)] + ['X', 'Y']
    return chrm_list

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_table', required=True, type=argparse.FileType('r'),
                        help='Tab delimited file of Sequenza segs and associated sample name.')
    parser.add_argument('-m', '--marker_file', required=True,
                        type=argparse.FileType('w'), help='Specify marker file output.')
    parser.add_argument('-o', '--out_seg', required=True,
                        type=argparse.FileType('w'), help='Specify GISTIC seg file output.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
