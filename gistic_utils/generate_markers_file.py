import argparse


def main():
    args = parse_args()
    window_size = args.window_size
    marker_str = "WINDOW_{1}_{0}\t{1}\t{2}"
    print "Marker Name\tChromosome\tMarker Position"
    for chrm in chrm_order:
        chrm_size = chrm_sizes[chrm]
        i = 0
        marker_pos = 1
        while marker_pos < chrm_size:
            print marker_str.format(i, chrm, marker_pos)
            marker_pos += window_size
            i += 1
    return


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--window_size", "-w", type=int, default=100000,
                        help='Size of window marker.')
    args = parser.parse_args()
    return args

chrm_sizes = {"1": 249250621,
              "2": 243199373,
              "3": 198022430,
              "4": 191154276,
              "5": 180915260,
              "6": 171115067,
              "7": 159138663,
              "8": 146364022,
              "9": 141213431,
              "10": 135534747,
              "11": 135006516,
              "12": 133851895,
              "13": 115169878,
              "14": 107349540,
              "15": 102531392,
              "16": 90354753,
              "17": 81195210,
              "18": 78077248,
              "19": 59128983,
              "20": 63025520,
              "21": 48129895,
              "22": 51304566,
              "X": 155270560,
              "Y": 59373566}

chrm_order = [str(x) for x in range(1, 23)] + ["X", "Y"]

if __name__ == '__main__':
    main()
