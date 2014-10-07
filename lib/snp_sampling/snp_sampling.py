#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import argparse
import textwrap
import random

def version():
    v1 = """
##########################################################################################
 Jinliang Yang
 updated: Oct.24.2013
 --------------------------------
 A step of XP-GWAS project
##########################################################################################
"""
    return v1

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
    parser.add_argument('-o', '--output', help='output files', type=str)
    
    parser.add_argument('--perc', help='100 percent of data to sample', default=50, type=int)
    parser.add_argument('--refcol', help='col of ref cumsum', default=5, type=int)
    parser.add_argument('--altcol', help='col of alt cumsum', default=6, type=int)
    parser.add_argument('--idcol', help='col of snp id', default=0, type=int)
    
    return parser
    #parser = get_parser()
    #parser.print_help()

def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

##########################################################################################

def read2dicts(file_name, refcol=7, altcol=8, idcol=10):
    
    start = 1
    tot = {}
    infile = open(file_name, 'r')
    next(infile) # skip the header line
    for line in infile:
        line = line.strip() # remove the newline and white spaces
        tokens = line.split()
        for idx in range(start, int(float(tokens[refcol]))+1):
            ref[idx] = tokens[idcol]
        start = int(float(tokens[refcol])) +1
        snphash[tokens[idcol]] = [0,0]
    infile.close()
    
    tot['ref'] = start-1
    ############################################    
    infile = open(file_name, 'r')
    next(infile)
    for line in infile:
        line = line.strip()
        tokens = line.split()
        for idx in range(start, int(float(tokens[altcol]))+1):
            alt[idx] = tokens[idcol]
        start = int(float(tokens[altcol])) +1 
    infile.close()
    tot['alt'] = start -1    

    return tot
#test: test = read2lofd(f)


        
def sampling():
    num = round(tot['alt']*args['perc']/100, 0)
    samplenum = random.sample(xrange(tot['alt']), int(num))
    for number in samplenum:
        if number <= tot['ref'] and number > 0:
            snpid = ref[number]
            snphash[snpid][0] += 1
        elif number > tot['ref']:
            snpid = alt[number]
            snphash[snpid][1] += 1

def writeFile():
	
    outfile = open(args['output'], 'w')
    outfile.write('\t'.join(['snpid', 'ref_count', 'alt_count']) + '\n')
	
    for k, v in snphash.items():
        outfile.write('\t'.join([k, str(v[0]), str(v[1])]) + '\n')
    outfile.close()            
##########################################################################################    
# main
parser = get_parser()
args = vars(parser.parse_args())

if args['input'] is None:
    print(version())
else:
    os.chdir(args['path'])


ref = {}
alt = {}
snphash = {}    
    
print('---loading the file: [', args['input'], '] ...')    
tot = read2dicts(args['input'], refcol=args['refcol'], altcol=args['altcol'], idcol=args['idcol'])

print('---sampling: [', args['perc'], '%] ...')
sampling()

print('---writing: [', args['output'], '] ...')
writeFile()