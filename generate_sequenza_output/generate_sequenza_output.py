#Written by Jasleen Grewal (jgrewal@bcgsc.ca)
#Date Created: 27th August 2015
#Last Modified: 14th September 2015

import subprocess
import commands
import os
import re
import argparse
import sys
import errno
import ntpath
from subprocess import call
import datetime
import time

def make_output_dir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def doesnt_exist_output(myfile):
    if os.path.isfile(myfile):
        return 0
    else:
        return 1

def main():
    """ Run sequenza_process.py
    """
    #Specify command line arguments
    parser = argparse.ArgumentParser(description="Generate CNV calls for input bam using Sequenza")
    parser.add_argument('-r', '--ref_fa', nargs=1, type=str, required=False, default = ["/reference/genomes/human/hg19a/genome/GRCh37-lite.fa"], help="Reference genome.fa file for mpileup step")
    parser.add_argument('-n', '--input_normal_bam', nargs=1, type=str, required=True, help="Input normal bam")
    parser.add_argument('-t', '--input_tumour_bam', nargs=1, type=str, required=True, help="Input tumour bam")
    parser.add_argument('-odir', '--output_dir', nargs=1, type=str, required=True, help="Output directory address")
    parser.add_argument('-su', '--seq_utils', nargs=1, type=str, required=False, default = ["/genesis/extscratch/morinlab/software/R_libs/sequenza/exec/sequenza-utils.py"],help="Sequenza utils python script location")
    parser.add_argument('-sa', '--seq_analyze', nargs=1, type=str, required=False, default=["/genesis/extscratch/morinlab/software/sequenza/run-sequenza/sequenza_analysis.R"], help="Sequenza analysis script")
    parser.add_argument('-gr', '--gc_ref', nargs=1, type=str, required=False, default = ["/genesis/extscratch/morinlab/shared/slin/GRCh37.gc50Base.txt.gz"], help="GC reference file for normalization")
    parser.add_argument('-s', '--sample_id', nargs=1, type=str, required=True, help="Name of sample (patient)")
    parser.add_argument('-g', '--gender', nargs=1, type=str, required=True, help="Gender of sample (patient). m or f")
    parser.add_argument('-p', '--parallel', default=False, help="Generate mpileups for each chromosome - reduce computational time",action="store_true")

    #Parse command line arguments
    args = parser.parse_args()
    input_fa = args.ref_fa[0]
    input_tumour = args.input_tumour_bam[0]
    input_normal = args.input_normal_bam[0]
    odir = os.path.abspath(args.output_dir[0])
    seq_utils = args.seq_utils[0]
    gc_ref = args.gc_ref[0]
    seq_analyze = args.seq_analyze[0]
    patient = args.sample_id[0]
    patient_gender = args.gender[0]
    is_parallel = 'False'#args.parallel
    #Create output dir(s)

    make_output_dir(odir+"/mpileups/"+patient)
    make_output_dir(odir+"/seqz/"+patient)
    make_output_dir(odir+"/seq_output/"+patient)
    final_odir =odir+"/seq_output/"+patient
    #Process input
    #mpileup object (1)
    my_command = []
    out_normal_mp = "{0}/mpileups/{1}/{2}.mpileup.gz".format(odir,patient,ntpath.basename(input_normal).replace(".bam",""))
    out_tumour_mp = "{0}/mpileups/{1}/{2}.mpileup.gz".format(odir,patient,ntpath.basename(input_tumour).replace(".bam",""))
    mp_normal = "samtools mpileup -f {0} -Q 20 {1} | gzip > {2}".format(input_fa,os.path.abspath(input_normal),out_normal_mp)
    mp_tumour = "samtools mpileup -f {0} -Q 20 {1} | gzip > {2}".format(input_fa,os.path.abspath(input_tumour),out_tumour_mp)

    if((doesnt_exist_output(out_normal_mp))):
        print "Will generate mpileup for .... " + input_normal.split("/",5)[-1]
        my_command.append(mp_normal)
    if((doesnt_exist_output(out_tumour_mp))):
        print "Will generate mpileup for .... " + input_tumour.split("/",5)[-1]
        my_command.append(mp_tumour)

    # seqz object (2)
    out_patient_seq ="{0}/seqz/{1}/{2}_tumour_seqz.gz".format(odir,patient,ntpath.basename(input_normal).replace(".bam",""))
    out_patient_seq_bin ="{0}/seqz/{1}/{2}_tumour_binned.seqz.gz".format(odir,patient,ntpath.basename(input_normal).replace(".bam",""))

    temp_seqz = "python {0} pileup2seqz -gc {1} -n {2} -t {3} | gzip > {4} ".format(seq_utils,gc_ref,out_normal_mp,out_tumour_mp,out_patient_seq)
    temp_seqz_bin = "python {0} seqz-binning -w 50 -s {1} | gzip > {2} ".format(seq_utils,out_patient_seq,out_patient_seq_bin)

    if((doesnt_exist_output(out_patient_seq))):
        print "Will generate joint normal_tumour binned (-w 50) seqz object for .... " + out_normal_mp.split("/",5)[-1] + " and " + out_tumour_mp.split("/",5)[-1]
        my_command.append(temp_seqz)

    if((doesnt_exist_output(out_patient_seq_bin))):
        print "Will generate joint normal_tumour binned (-w 50) seqz object for .... " + out_normal_mp.split("/",5)[-1] + " and " + out_tumour_mp.split("/",5)[-1]
        my_command.append(temp_seqz_bin)

    # R analysis object (3)
    temp_r = "Rscript {0} --input {1} --output {2} --sample_id {3} --gender {4} --ploidy_limit 2 --min_reads 20 --max_cn 5".format(seq_analyze,out_patient_seq_bin,final_odir,patient,patient_gender)
    print "Will analyze seqz object and generating R analysis output at .... " + final_odir
    my_command.append(temp_r)

    my_command = "  && ".join(my_command)
    #print my_command
    p = subprocess.Popen(my_command,stdout=subprocess.PIPE,shell=True) #os.system(my_command)
    log_output = open(os.path.join(odir, "log_process.log"), 'w+')
    log_output.write("\n------------\n {0}\n--------\n processed sample {1} with command {2}".format(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), patient,my_command))
    out, err = p.communicate()
    log_stderr = open(os.path.join(odir, "log_stderr.log"), 'w+')
    log_stderr.write(str(err))
    log_stdout = open(os.path.join(odir, "log_stdout.log"), 'w+')
    log_stdout.write(str(out))
    #os.system(my_command)
    #log_output = open(odir+"log_process.log",'w+'); log_output.write("\n------------\n {0}\n--------\n processed sample {1} with command {2}".format(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), patient,my_command))

if __name__ == '__main__':
    main()
