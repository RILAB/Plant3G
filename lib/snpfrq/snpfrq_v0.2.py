#!/usr/bin/env python

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
    snpfrq version 1.0
    Jinliang Yang
    updated: July.3.2014
    --------------------------------

    compute SNP frq and loci missing rate from DSF
    ##########################################################################################
    """
    return v1



def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

##########################################################################################

def readfile_and_process(infile_name, outfile_name):

    with open(infile_name, 'r') as infile:
    
        line1 = infile.readline()
        line1array = line1.split()
        if(line1array[0] != "snpid"):
            warning("First col should be snpid!")
        line2 = infile.readline()
        line2array = line2.split()
        if(len(str(line2array[start])) != 1):
            warning("SNP should be coded with one latter, like:'A T C G or - +'!")

    with open(infile_name, 'r') as infile, open(outfile_name, "w") as outfile:

        line1 = infile.readline()
        line1data = line1.split()
        outfile.write("\t".join([line1data[0],"major","minor","MAF","missing"]) + "\n")
        for line in infile:
            tokens = line.split()
            out = get_loci_info(tokens[start:end])

            ### print out the results
            if out:
                outfile.write("\t".join([tokens[0], out['major'], out['minor'], str(out['maf']), str(out['missing'])]) + "\n")

#test: test = read2lofd(f)

def write_prob_snp():
    """
    Write the problem snp to the input.prob file:
    """
    output = args['output'].split('.')[0]
    output = ".".join([output, "mul"])
    outfile = open(output, 'w')
    for prob_snp in prob:
        outfile.write('\t'.join(prob_snp) + '\n')
    outfile.close()



def get_loci_info(tokens):
    
    set0 = set(tokens)
    if 'N' in set0:
        set0.remove('N')
    snpset = list(set0)
        
    info = {}
    if len(snpset) > 2:
        prob.append(tokens)
    elif len(snpset) == 2:
        c1 = tokens.count(snpset[0])
        c2 = tokens.count(snpset[1])
   
        if c1 >= c2:
            info['major'] = snpset[0]
            info['minor'] = snpset[1]
            info['maf'] = round(c2/(c1+c2),3)
        else:
            info['major'] = snpset[1]
            info['minor'] = snpset[0]
            info['maf'] = round(c1/(c1+c2),3)
        info['missing'] = round((len(tokens) - c1 - c2)/len(tokens),3)

    return info

##########################################################################################
#get_loci_info(y[4:31])
    
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
    parser.add_argument('-s','--start', help='start cols (1-based) of the genotype', type=int)
    parser.add_argument('-e','--end', help='end cols (1-based) of the genotype', type=int)

    parser.add_argument('-o', '--output', help='output files, like chr1_merged', type=str)

    return parser
    #parser = get_parser()
    #parser.print_help()

if __name__ == '__main__':
    parser = get_parser()
    args = vars(parser.parse_args())

    if args['input'] is not None:
        print(version())
    if args['path'] is not None:
        os.chdir(args['path'])

    ##### read in the input file ######
    st = timeit.default_timer()
    prob=[]
    print("Reading and writing SNP info ...")

    start = args['start'] -1
    end = args['end']
    readfile_and_process(args['input'], args['output'])

    print("Writing SNPs with multiple alleles ...")
    write_prob_snp()

    et = timeit.default_timer()

    print("[ ", "%.0f" % ((et - st)/60), " ] minutes of run time!")
    print("Job finished!")