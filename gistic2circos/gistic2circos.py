import argparse
import re

def main():
    args = parse_args()
    all_lesion_file = args.all_lesion_file

    header_found = False

    peak_dict = {}

    chrm_order = [ str(c) for c in range(1, 23) ]

    for peak in all_lesion_file:

        if not header_found:
            header_found = True
            continue

        fields = peak.split('\t')

        if "CN values" in fields[0]:
            continue

        is_amp = False
        cnv_state = None
        chrm = None
        chrm_start = None
        chrm_end = None

        if "Amplification" in fields[0]:
            is_amp = True

        p = re.compile("chr(\d+):(\d+)-(\d+)\(.+\)")
        m = re.search(p, fields[2])

        chrm = m.group(1)
        chrm_start = m.group(2)
        chrm_end = m.group(3)

        if is_amp:
            cnv_state = 'amp'
        else:
            cnv_state = 'homd'

        if chrm not in peak_dict:
            peak_dict[chrm] = []

        peak_dict[chrm].append([chrm_start, chrm_end, cnv_state])

        #print "{}\t{}\t{}\t{}".format(chrm, chrm_start, chrm_end, cnv_state)

    for chrm in chrm_order:
        if chrm in peak_dict:
            for peak in peak_dict[chrm]:
                print "{}\t{}\t{}\t{}".format(chrm, peak[0], peak[1], peak[2])

    return

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('all_lesion_file', type=argparse.FileType('r'),
                        help="all lesions result files from GISTIC output. E.g., 'all_lesions.conf_90.txt'")

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    main()
