import pandas
import argparse
import numpy
import os

def main():

    args = parse_args()

    input_files = args.input
    output_filenames = args.output_filenames
    is_genefilter = args.genebased
    header_indices = args.header_row
    outdir = args.output_dir[0]

    if len(input_files) < 2:
        print "ERROR: Specify at least two MAF files."
        return

    if len(input_files) != len(output_filenames):
        print "ERROR: Number of output filenames do not match the number of input MAFs"
        return

    if not len(header_indices):
        header_indices = [0] * len(input_files)
    elif len(input_files) != len(header_indices):
        print "ERROR: Number of header indices do not match the number of input MAFs"
        return

    if not os.path.exists(outdir):
        print "Making directory: " + outdir
        os.makedirs(outdir)

    parsed_mafs = parse_mafs_into_array(input_files,header_indices)

    while len(parsed_mafs) > 0:
        maf = parsed_mafs.pop(0)
        for i in range(len(parsed_mafs)):
            parsed_mafs[i] = onlyInLeft(parsed_mafs[i],maf,is_genefilter)

        output_filename = output_filenames.pop(0)
        maf.to_csv((outdir + "/" + output_filename),sep="\t",index=False)
        print "Variants collapsed to " + outdir + "/" + output_filename

    return

def parse_mafs_into_array(input_files,header_indices):
    parsed_mafs = []
    i = 0
    for input_file in input_files:
        parsed_maf = pandas.DataFrame.from_csv(input_file,sep="\t",index_col=None,header=header_indices[i])
        parsed_maf = updatecoltypes(parsed_maf)
        parsed_mafs.append(parsed_maf)
        i += 1

    return parsed_mafs

def onlyInLeft(table1,table2,genefilter_flag):
    if (genefilter_flag):
        return table1[~((table1.Hugo_Symbol.isin(table2.Hugo_Symbol)))]
    else:
        return table1[~((table1.Chromosome.isin(table2.Chromosome)) & (table1.Start_Position.isin(table2.Start_Position)) & (table1.End_Position.isin(table2.End_Position)))]

def updatecoltypes(mytable):
    temptable=makecolint(mytable,["Start_Position","End_Position"])
    temptable=makecolstr(temptable,["Chromosome"])
    return temptable

def makecolint(table1,listofcols):
    table1[listofcols] = table1[listofcols].astype(numpy.int64)
    return table1

def makecolstr(table1,listofcols):
    table1[listofcols] = table1[listofcols].astype(str)
    return table1

def parse_args():
    parser = argparse.ArgumentParser(description="Collapses shared variants within MAF files")

    parser.add_argument('-i','--input',nargs='+',type=argparse.FileType('r'),required=True,help='Specify two or more MAFs.')
    parser.add_argument('-o','--output_filenames',nargs='+',type=str,required=True,help='Specify output filenames for unique MAFs. Number of prefixes must equal number of input MAFs. Order of filenames should correspond to order of input MAFs.')
    parser.add_argument('-g','--genebased',action='store_true',default=False,help='Instead of filtering by exact mutation location, filter by gene name.')
    parser.add_argument('-r','--header_row',nargs='+',type=int,help='If your input files have different header indices, pass them here. Number of header indices must equal number of input MAFs. Defaults to first line in all MAFs.')
    parser.add_argument('-d','--output_dir',nargs=1,type=str,required=True,help='Specify output directory for unique MAFs.')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    main()
