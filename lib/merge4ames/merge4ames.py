#!/usr/bin/env python

from __future__ import print_function
from collections import Counter
import sys
import argparse
import textwrap
import os

def version():
    v1 = """
    ##########################################################################################
    merge4ames version 1.0
    Jinliang Yang
    updated: 2.11.2014
    --------------------------------

    Merge replicated samples of the Ames panel GBS data to produce consensus SNP callings.
    When comparing two genotypes, only replace the genotype in the original list
    if it is missing and the second genotype is not or make it missing if the two genotypes
    do not agree!
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
    parser.add_argument('--cutoff', help='cutoff of the consensus', type=float)
    parser.add_argument('--idx', help='path of the idx file', type=str)
    parser.add_argument('--infile', help='path of the input chr data', type=str)
    parser.add_argument('-f', '--infiles', help='a list of merging files', type=str)
    parser.add_argument('-o', '--output', help='output files, like chr1_merged', type=str)
    parser.add_argument('-m','--mode',
                        help=''' mode for output format:
                        [--mode=1 default],
                        0, Other not determined;
                        1, GenSel format;
                        ''',
                        choices=[0,1], default=1, type=int)
    return parser

### use print_function
def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

##########################################################################################
def read_idx(ames_idx_file):
    temp = []
    infile = open(ames_idx_file, 'r')
    next(infile)
    for line in infile:
        line = line.strip()
        token = line.split('\t')
        temp.append(token)
    infile.close()
    return temp # return a list of list ['PI405705', '16 17 18 19 20 21', '6']

#The GBS genotype files in HapMap format (*.hmp.txt.gz) can be opened directly in TASSEL in
#their compressed (*.gz) state. To save space, we use single letters to encode unphased,
#diploid genotypes, where:
#A = A/A
#C = C/C
#G = G/G
#T = T/T
#K = G/T
#M = A/C
#R = A/G
#S = C/G
#W = A/T
#Y = C/T
#N = missing
def read_snp(snpfile):
    temsnp = []
    infile = open(snpfile, 'r')
    fline = infile.readline()
    nms = fline.split()

    for line in infile:
        line = line.strip()
        token = line.split()
        for i in xrange(11, len(token)):
            if token[i] == 'A':
                token[i] = 'A A'
                count['ATCG'] += 1
            elif token[i] == 'C':
                token[i] = 'C C'
                count['ATCG'] += 1
            elif token[i] == 'G':
                token[i] = 'G G'
                count['ATCG'] += 1
            elif token[i] == 'T':
                token[i] = 'T T'
                count['ATCG'] += 1
            elif token[i] == 'K':
                token[i] = 'G T'
                count['Y'] += 1
            elif token[i] == 'M':
                token[i] = 'A C'
                count['Y'] += 1
            elif token[i] == 'R':
                token[i] = 'A G'
                count['Y'] += 1
            elif token[i] == 'S':
                token[i] = 'C G'
                count['Y'] += 1
            elif token[i] == 'W':
                token[i] = 'A T'
                count['Y'] += 1
            elif token[i] == 'Y':
                token[i] = 'C T'
                count['Y'] += 1
            elif token[i] == 'N':
                token[i] = 'N N'
                count['N'] += 1
            else:
                warning("Non reconginzed alleles found!", token[i])
        temsnp.append(token)
    infile.close()
    return(temsnp)

def get_maf_missing(asnp=snp[0][11:]):
    snptable = Counter(asnp)
    dict = {}
    if 'N N' in snptable:
        dict['missing'] = float(snptable['N N'])/(len(asnp) - 11)
        del(snptable['N N'])
    else:
        dict['missing'] = 0
     
    if len(snptable) == 1:
        dict['maf'] = 0
        for key in snptable:
            dict['major'] = key.split()[0]
            dict['minor'] = dict['major']
    elif len(snptable) == 2 or len(snptable) == 3:
        mysnp = {} #assign a hash to A or T SNP type
        for key in snptable:
            temsnp = key.split()
            if temsnp[0] == temsnp[1]:
                if temsnp[0] in mysnp:
                    mysnp[temsnp[0]] = mysnp[temsnp[0]] + 2*snptable[key]
                else:
                    mysnp[temsnp[0]] = 2*snptable[key]
            elif temsnp[0] != temsnp[1]:
                if temsnp[0] in mysnp:
                    mysnp[temsnp[0]] = mysnp[temsnp[0]] + snptable[key]
                else:
                    mysnp[temsnp[0]] = snptable[key]
                if temsnp[1] in mysnp:
                    mysnp[temsnp[1]] = mysnp[temsnp[1]] + snptable[key]
                else:
                    mysnp[temsnp[1]] = snptable[key]
        #### end of the for loop
        if len(mysnp) == 2:
            ### sort: A T
            sortedsnp = sorted(mysnp.values())
            dict['maf'] = float(sortedsnp[0])/(sortedsnp[0] + sortedsnp[1])
            sortedkey = sorted(mysnp, key=mysnp.__getitem__)
            dict['minor'] = sortedkey[0]
            dict['major'] = sortedkey[1]
        else:
            dict['maf'] = -9 ### multiple alleles
            for key in mysnp:
                dict['minor'] = ' '.join(key)
                dict['major'] = ' '.join(key)
    else:
        dict['maf'] = -8 ### multiple alleles
        for key in snptable:
            dict['minor'] = ' '.join(key)
            dict['major'] = ' '.join(key)
    return dict

def consensus_snp_call(idxob=amesidx[5][1], onesnp=snp[0], CUTOFF=0.6):

    idxlist = idxob.split()
    snpcalls =[]
    for idx in idxlist:
        snpcalls.append(onesnp[int(idx) -1])
    snptable2 = Counter(snpcalls)
    if 'N N' in snptable2:
        del(snptable2['N N'])

    if len(snptable2) == 0:
        mysnp = ['N N']
    elif len(snptable2) == 1:
        mysnp = snptable2.keys()
    elif len(snptable2) >= 2:
        sortedsnp = sorted(snptable2, key=snptable2.__getitem__)
        sortedval = sorted(snptable2.values())
        frq = float(sortedval[-1])/sum(sortedval)
        if frq >= CUTOFF:
            mysnp = [sortedsnp[-1]]
        else:
            mysnp = ['N N']
    else:
        mysnp = 'ERROR'
        print('Error detected for SNP:', onesnp[0])
    return mysnp

#snp[3][15:20]
#consensus_snp_call(idxob=amesidx[4][1], onesnp=snp[3])

def one_snp_merge(amesidx=amesidx, onesnp=[]):

    """

    :rtype : resout
    """
    snpinfo = get_maf_missing(asnp= onesnp[11:])
    if snpinfo['maf'] == 0:
        # SNP has no polymorphism
        resout['noploy'].append(onesnp[0])
    elif snpinfo['maf'] < 0:
        # SNP have multiple alleles
        resout['multi'].append(onesnp[0])
    else:
        # merge one SNP for each accession
        temsnp = [onesnp[0]]
        for plant in amesidx:
            if plant[2] == 1:
                temsnp.extend(onesnp[int(plant[1]) - 1])
            # merge duplicated calls with consensus mode
            elif plant[2] > 0:
                snpcall = consensus_snp_call(idxob=plant[1], onesnp=onesnp, CUTOFF=0.6)
                temsnp.extend(snpcall)
            else:
                print("Error type 1: amesidx error", amesidx[i])
        aftermaf = get_maf_missing(asnp=temsnp[1:])
        resout['snp'] = [temsnp[0]] + [aftermaf['minor']] + [aftermaf['maf']] + [aftermaf['missing']] + temsnp[1:]

def merge_all():
    # results:
    count = {'ATCG':0, 'Y':0, 'N':0}
    resout = {'noploy':[], 'multi':[], 'snp':[]}

    print("---- reading index file:", os.path.basename(args['idxfile']));
    idx = read_idx(args['idxfile'])

    snpall = []
    #### process all file together for GenSel:
    if args['batch'] == 1:
        infile = open(args['infiles'], 'r')
        for line in infile:
            chr = read_snp(snpfile=line)
            print("- reading: [", len(chr), "] SNPs from: [", os.path.basename(line), "]")
            snpall.append(chr)

    #### process one chr at a time
    if args['batch'] == 0:
        chr = read_snp(snpfile=args['infile'])
        print("- reading: [", len(chr), "] SNPs from: [", os.path.basename(args['infile']), "]")
        snpall = chr

    snpout = []
    for snpi in xrange(snpall):
        one_snp_merge(amesidx=idx, onesnp=snpi)
    mode = args['mode']
    if mode == 1:
        for mysnp in xrange(resout['snp']):
            for i in range(4,len(mysnp)):
                twocall = mysnp[i].split()
                if twocall[0] != twocall[1]:
                    mysnp[i] = 0
                elif twocall[0] == "N":
                    mysnp[i] = 0
                elif twocall[0] == mysnp[]




def recoding_output(snpline):

    if mode == 1: #GenSel

        #### SNP transpose
        for snpi in snpall:
            temsnp = snpi.split()[0]
            if temsnp == minor: # -10
                snpout.append(-10)
            elif temsnp != minor: #10
                snpout.append(10)
            elif temsnp == "N":
                snpout.append(0)






            for snpj in snpi:

            gensel[1]

    outfile = open(args['output'], "w")
    outfile.write("\t".join(header) + )
    for i in res:
        outfile.write("\t".join(snp) + '\n')
    outfile.close()




def writeLog():

    outfile1 = open(args['output'] + '.noploy', 'w')
    outfile1.write('\n'.join(resout['noploy']))
    outfile1.close()
    
    outfile2 = open(args['output'] + '.multi', 'w')
    outfile2.write('\n'.join(resout['multi']))
    outfile2.close()
    
    outfile3 = open(args['output'] + '.log', 'w')
    outfile3.write()


##########################################################################################


def writeLog():
	"""
	Write a log file containing the initial SNP counts, the number of merged SNPs, and the final SNP count.

	Input: A string containing the path to write the log to.
	"""

	outfile = open(args['output'] + '.log', 'w')

	outfile.write('---Initial SNPs:---\n')
	outfile.write('Hapmap1: ' + str(len(hapmap1)) + '\nHapmap2: ' + str(len(hapmap2)) \
	 + '\nRNA-seq: ' + str(len(rnaseq)) + '\n')

	total = len(hapmap1) + len(hapmap2) + len(rnaseq)
	outfile.write('Total: ' + str(total) + '\n')

	merged_snps = total - len(res)
	outfile.write('---Merged SNPs: ' + str(merged_snps) + '\n')
	outfile.write('---Final SNPs: ' + str(len(res)) + '\n')
	outfile.write('\n')

	outfile.write('---Dataset unique loci---\n')
	outfile.write('-HapMap1 only: ' + str(len(unq['hmp1'])) + '\n')
	outfile.write('-HapMap2 only: ' + str(len(unq['hmp2'])) + '\n')
	outfile.write('-RNA-seq only: ' + str(len(unq['rnaseq'])) + '\n')

	outfile.write('---Dataset matched loci and alleles ---\n')
	outfile.write('-HapMap1 and HapMap2: ' + str(len(lmat['hmp12'])) + '\t' + str(amat['hmp12']) + '\n')
	outfile.write('-HapMap1 and RNA-seq: ' + str(len(lmat['hmp1rna'])) + '\t' + str(amat['hmp1rna']) + '\n')
	outfile.write('-HapMap2 and RNA-seq: ' + str(len(lmat['hmp2rna'])) + '\t' + str(amat['hmp2rna']) + '\n')
	outfile.write('-HapMap1 and HapMap2 and RNA-seq: ' + str(len(lmat['hmp12rna'])) + '\t' + str(amat['hmp12rna']) + '\n')

	outfile.write('---Dataset unmatched loci and alleles ---\n')
	outfile.write('-HapMap1 and HapMap2: ' + str(len(lunmat['hmp12'])) + '\t' + str(aunmat['hmp12']) + '\n')
	outfile.write('-HapMap1 and RNA-seq: ' + str(len(lunmat['hmp1rna'])) + '\t' + str(aunmat['hmp1rna']) + '\n')
	outfile.write('-HapMap2 and RNA-seq: ' + str(len(lunmat['hmp2rna'])) + '\t' + str(aunmat['hmp2rna']) + '\n')
	outfile.write('-HapMap1 and HapMap2 and RNA-seq: ' + str(len(lunmat['hmp12rna'])) + '\t' + str(aunmat['hmp12rna']) + '\n')

	outfile.close()

########################################################################################################################
'''
main
'''
parser = get_parser()
args = vars(parser.parse_args())

if args['infiles'] or args['hmp1'] is None:
    print(version())

if args['path'] is not None:
    os.chdir(args['path'])

header = ['rs', 'alleles', 'chr', 'pos', 'B73', 'Z001', \
'Z002', 'Z003', 'Z004', 'Z005', 'Z006', 'Z007', 'Z008', 'Z009', \
'Z010', 'Z011', 'Z012', 'Z013', 'Z014', 'Z015', 'Z016', 'Z017', \
'Z018', 'Z019', 'Z020', 'Z021', 'Z022', 'Z023', 'Z024', 'Z025', \
'Z026', 'id' ]
hapmap1 = []
hapmap2 = []
rnaseq = []

#### read in the SNP data ####
if args['infiles'] is not None:
    input = read_batch()
else:
    input = read_chr()

#### output statistics ####
res = []
unq = {'hmp1':[], 'hmp2':[], 'rnaseq':[]}
lmat = {'hmp12':[], 'hmp1rna':[], 'hmp2rna':[], 'hmp12rna':[]}
lunmat = {'hmp12':[], 'hmp1rna':[], 'hmp2rna':[], 'hmp12rna':[]}
amat = {'hmp12':0, 'hmp1rna':0, 'hmp2rna':0, 'hmp12rna':0}
aunmat = {'hmp12':0, 'hmp1rna':0, 'hmp2rna':0, 'hmp12rna':0}


print('Merging the SNP data...')
merge_snp()

print('Writing merged data to file...')
writeMergedData()

print('Writing log file...')
writeLog()
