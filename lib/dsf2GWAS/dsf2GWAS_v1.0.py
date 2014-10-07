#!/usr/bin/env python
#__author__ = 'yangjl'

from __future__ import print_function
from __future__ import division
import sys
import argparse
import textwrap
import timeit
import os




def version():
    v1 = """
    ##########################################################################################
    version 1.0
    Jinliang Yang
    updated: July 9th, 2014
    passed R code test
    --------------------------------

    transform density SNP format (dsf) to GWAS (gensel4.2 now) and other formats
    ##########################################################################################
    """
    return v1

def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

##########################################################################################
def checkfile(infile_name):

    with open(infile_name, 'r') as infile:
        ### check the first line
        line1 = infile.readline()
        line1a = line1.split()
        if(line1a[0] != "snpid" or line1a[1] != "major" or line1a[2] != "minor"):
            warning("snpid major and minor should be the first three columns in the header")
        else:
            print("input file header OK!")

        ### check the 2nd line
        line2 = infile.readline()
        line2array = line2.split()
        if(len(str(line2array[start])) != 1):
            warning("SNP should be coded with one latter, like:'A T C G or - +'!")
        else:
            print("input file coding format OK!")

###############################
def readfile(infile_name):

    with open(infile_name, 'r') as infile:
        tsnp = []
        line1 = infile.readline()
        line1array = line1.split()
        tsnp.append(line1array)

        for line in infile:
            tokens = line.split()
            if mode == 0:
                tsnp.append(recode2gensel(asnp=tokens))
            else:
                warning("only support mode=0 currently!")

        return tsnp

################################
def writefile(outfile_name, tsnp=[]):
    with open(outfile_name, "w") as outfile:
        if mode == 0:
            #http://stackoverflow.com/questions/6473679/python-list-of-lists-transpose-without-zipm-thing
            # [[j[i] for j in l] for i in range(len(l))]
            zipsnp = map(list, zip(*tsnp))

            outfile.write('\t'.join(zipsnp[0]) + '\n')
            for snps in zipsnp[start:end]:
                outfile.write('\t'.join(snps) + '\n')


#################################
def recode2gensel(asnp=[]):
    '''
    diallel imputation for GenSel: major=10, minor=-10, missing, heter=0
    '''
    mysnp = asnp[0:start]
    #print(len(asnp))
    # major = 10, minor = -10
    major = asnp[1]
    minor = asnp[2]

    for i in range(start, end):
        if asnp[i] == "N":
            mysnp.append("0")
        elif asnp[i] == major:
            mysnp.append("10")
        elif asnp[i] == minor:
            mysnp.append("-10")
        else:
            warning(asnp[i], "have multiple alleles when impute for gensel!")
    return mysnp



################################
def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(version())
        )

    # positional arguments:
    #parser.add_argument('query', metavar='QUERY', type=str, nargs='*', \
    #                    help='the question to answer')

    # optional arguments:
    parser.add_argument('-p', '--path', help='the path of the input files', \
                        nargs='?', default=os.getcwd())
    parser.add_argument('-i','--input', help='input file', type=str)
    parser.add_argument('-o', '--output', help='output files, like chr1_merged', type=str)
    parser.add_argument('-s','--start', help='start position (1-based) of the genotype', type=int)
    parser.add_argument('-e', '--end', help='end position (1-based) of the genotype', type=int)
    parser.add_argument('-m','--mode',
                        help='''[default --mode=0];
                        0, for GenSel [currently support];
                        1, for PLINK;
                        2, for SNPTEST;
                        ''',
                        choices=[0, 1,2], default=0, type=int)
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = vars(parser.parse_args())

    if args['input'] is not None:
        print(version())
    if args['path'] is not None:
        os.chdir(args['path'])

    ##### read in the pedigree file ######
    st = timeit.default_timer()

    ##### read in the density snp file ######
    start = args['start']-1
    end = args['end']
    mode = args['mode']

    checkfile(infile_name = args['input'])

    print("Reading density SNP format (dsf) file ...")
    dsnp = readfile(infile_name = args['input'])
    print("[ ", len(dsnp), " ] SNPs loaded from [ ", end-start, " ] founders!")

    ##### start imputation ######
    print("Start recoding using mode [ %d ]" % mode)
    writefile(args['output'], tsnp=dsnp)
    print("File output to: [ %s ]" % args['output'])

    et = timeit.default_timer()
    print("[ %.0f ] minutes of run time!" % ((et - st)/60))
    print("Job Done!")




