#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
import sys
import argparse
import textwrap
import timeit
import os



def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

##########################################################################################
def read_ped(pedfile):
    temp = []
    infile = open(pedfile, 'r')
    next(infile) # skip the header line
    for line in infile:
        tokens = line.split()
        temp.append(tokens)
    infile.close()
    return temp

##########################################################################################
def checkFile(infile_name):

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
        for snpi in range(start, end):
            if(len(str(line2array[snpi])) != 1):
                warning("SNP should be coded with one latter, like:'A T C G or - +'!")
        else:
            print("input file coding format OK!")

def readLofD(infile_name):
    temp = []
    with open(infile_name, 'r') as infile:
        next(infile) # skip the header line
        for line in infile:
            tokens = line.split()
            key = tokens[0]
            temp.append({key:tokens})
        return temp

for i in range(0, len(first)):
            temp1[first[i]] = tokens[i]
##########################################################################################
def impute_raw(loci=[], ped=[]):
    '''
    raw type
    '''
    mysnp = [loci['chr'], loci['id'], loci['pos']]
    for aped in ped:
        psnp1 = loci[aped[0]]
        psnp2 = loci[aped[1]]
        # using original SNP coding "ATCG"
        major = loci['major']
        minor = loci['minor']
        
        # probility of AA, AB and BB        
        if psnp1 == major and psnp2 == major:
            mysnp.append('\t'.join([major, major]))
        elif psnp1 == major and psnp2 == minor:
            mysnp.append('\t'.join([major, minor]))
        elif psnp1 == major and psnp2 == "N":
            mysnp.append('\t'.join([major, 'N']))
        elif psnp1 == minor and psnp2 == major:
            mysnp.append('\t'.join([major, minor]))
        elif psnp1 == minor and psnp2 == minor:
            mysnp.append('\t'.join([minor, minor]))
        elif psnp1 == minor and psnp2 == "N":
            mysnp.append('\t'.join([minor, 'N']))
        elif psnp1 == "N" and psnp2 == major:
            mysnp.append('\t'.join([major, 'N']))
        elif psnp1 == "N" and psnp2 == minor:
            mysnp.append('\t'.join([minor, 'N']))      
        elif psnp1 == "N" and psnp2 == "N":
            mysnp.append('\t'.join(['N', 'N']))
        else:
            warnings(loci['id'], "have gensel imputation error!!!")
    return mysnp   

    
def impute_gensel(aped=[], dsnp=[]):
    '''
    diallel imputation for GenSel: major=10, minor=-10, missing, heter=0
    '''
    mysnp = []
    for loci in dsnp:
        psnp1 = loci[aped[0]]
        psnp2 = loci[aped[1]]
        # major = 10, minor = -10
        major = loci['major']
        minor = loci['minor']
                
        if psnp1 == major and psnp2 == major:
            mysnp.append('10')
        elif psnp1 == major and psnp2 == minor:
            mysnp.append('0')
        elif psnp1 == major and psnp2 == "N":
            mysnp.append('5')
        elif psnp1 == minor and psnp2 == major:
            mysnp.append('0')
        elif psnp1 == minor and psnp2 == minor:
            mysnp.append('-10')
        elif psnp1 == minor and psnp2 == "N":
            mysnp.append('-5')
        elif psnp1 == "N" and psnp2 == major:
            mysnp.append('5')
        elif psnp1 == "N" and psnp2 == minor:
            mysnp.append('-5')      
        elif psnp1 == "N" and psnp2 == "N":
            mysnp.append('0')
        else:
            warning(loci['id'], "have gensel imputation error!!!")
    return mysnp


    
def impute_plink(aped=[], dsnp=[]):
    '''
    diallel imputation for PLINK: major="MM", minor="GG", missing, heter="NN"
    '''
    mysnp = []
    for loci in dsnp:
        psnp1 = loci[aped[0]]
        psnp2 = loci[aped[1]]
        # using original SNP coding "ATCG"
        major = loci['major']
        minor = loci['minor']
                    
        if psnp1 == major and psnp2 == major:
            mysnp.append(' '.join([major, major]) )
        elif psnp1 == major and psnp2 == minor:
            mysnp.append(' '.join([major, minor]) )
        elif psnp1 == major and psnp2 == "N":
            mysnp.append('0 0')  
        elif psnp1 == minor and psnp2 == major:
            mysnp.append(' '.join([major, minor]) )
        elif psnp1 == minor and psnp2 == minor:
            mysnp.append(' '.join([minor, minor]) )
        elif psnp1 == minor and psnp2 == "N":
            mysnp.append('0 0') 
        elif psnp1 == "N" and psnp2 == major:
            mysnp.append('0 0')
        elif psnp1 == "N" and psnp2 == minor:
            mysnp.append('0 0')      
        elif psnp1 == "N" and psnp2 == "N":
            mysnp.append('0 0')
        else:
            warning(loci['id'],psnp1, psnp2, "have plink imputation error!!!")
    return mysnp
    
def impute_snptest(loci=[], ped=[]):
    '''
    IMPUTE Diallel for SNPTESTv2
    The next three numbers on the line should be the probabilities of the three 
    genotypes AA, AB and BB at the SNP for the first individual in the cohort
    '''
    mysnp = [loci['chr'], loci['id'], loci['pos'], loci['major'], loci['minor']]
    for aped in ped:
        psnp1 = loci[aped[0]]
        psnp2 = loci[aped[1]]
        # using original SNP coding "ATCG"
        major = loci['major']
        minor = loci['minor']
        
        # probility of AA, AB and BB        
        if psnp1 == major and psnp2 == major:
            mysnp.append('1 0 0')
        elif psnp1 == major and psnp2 == minor:
            mysnp.append('0 1 0')
        elif psnp1 == major and psnp2 == "N":
            mysnp.append('0.5 0.5 0')
        elif psnp1 == minor and psnp2 == major:
            mysnp.append('0 1 0')
        elif psnp1 == minor and psnp2 == minor:
            mysnp.append('0 0 1')
        elif psnp1 == minor and psnp2 == "N":
            mysnp.append('0 0.5 0.5')
        elif psnp1 == "N" and psnp2 == major:
            mysnp.append('0.5 0.5 0')
        elif psnp1 == "N" and psnp2 == minor:
            mysnp.append('0 0.5 0.5')      
        elif psnp1 == "N" and psnp2 == "N":
            mysnp.append('0 0 0')
        else:
            warnings(loci['id'], "have gensel imputation error!!!")
    return mysnp     
    
##########################################################################################            
def write_prob_snp():
    """
    Write the problem snp to the input.prob file:
    """
    output = args['dsnp'].split('.')[0]
    output = ".".join([output, "prob"])
    outfile = open(output, 'w')
    for prob_snp in prob:
        outfile.write('\t'.join(prob_snp) + '\n')

    outfile.close()

def impute_write(ped=[], dsnp=[], mode=0):
    
    outfile = open(output, 'w')
    if args['header'] == "yes":
        outfile.write('ID' + '\t')
        for loci in dsnp:
            outfile.write('\t'.join(loci['id']) +'\n')
    
    if mode == 0:
        for aped in ped:
            print("impute:", len(dsnp), " SNPs at mode=gensel for [ ", aped[2], " ]", "\r")
            imputesnp = impute_gensel(aped=aped, dsnp=dsnp)
            outfile.write(aped[2] + '\t' + '\t'.join(imputesnp) + '\n')
    if mode == 1:
        for aped in ped:
            print("impute:", len(dsnp), " SNPs at mode=plink for [ ", aped[2], " ]", "\r")
            imputesnp = impute_plink(aped=aped, dsnp=dsnp)
            outfile.write('\t'.join([aped[2],'1','0','0','1','0']) + "\t")
            outfile.write('\t'.join(imputesnp) + '\n')
    if mode == 2:
        print("impute", len(dsnp), " SNPs at mode=snptest")
        for loci in dsnp:
            imputesnp = impute_snptest(loci=loci, ped=ped)
            outfile.write('\t'.join(imputesnp) + '\n')
    if mode == 3:
        print("impute", len(dsnp), " SNPs at mode=raw")
        for loci in dsnp:
            imputesnp = impute_snptest(loci=loci, ped=ped)
            outfile.write('\t'.join(imputesnp) + '\n')
    outfile.close()
                
def write_pheno(mode=0):

    if mode == 0:
        outpheno = ".".join([output, "pheno"])
        outfile = open(outpheno, "w")
        outfile.write('\t'.join(["Genotype", "Pheno", "Fix", "$cof"]) + '\n')
        for aped in ped:
            outfile.write('\t'.join([aped[2], '-9', '1', '0']) + '\n')
        outfile.close()
        
    if mode == 1:
        outpheno = ".".join([output, "map"])
        outfile = open(outpheno, "w")
        for loci in dsnp:
            outfile.write('\t'.join([loci['chr'], loci['id'], '0', loci['pos'] ]) + '\n')
        outfile.close()
            
    if mode ==2:
        outpheno  =".".join([output, "sample"])
        outfile = open(outpheno, "w")
        outfile.write('\t'.join(['ID_1', 'ID_2', 'missing', 'cov1', 'cov2', 'pheno1', 'bin1' ]) +'\n')
        outfile.write('\t'.join(['0', '0', '0', 'D', 'C', 'P', 'B' ]) +'\n')
        for aped in ped:
            outfile.write("\t".join([aped[2], '1', '0', '1', '2', '-9', '0']) + '\n')
    return outpheno    
          
    
def version():
    v = """
################################################################################
 impute4diallel version 0.2
 Jinliang Yang
 updated: July.12.2013
 changes: input changed to sdf
 --------------------------------

 SNP Imputation for diallel!
################################################################################
"""
    return v

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
    parser.add_argument('-d','--ped', help='pedigree of the F1 (header: nm, p1, p2)', type=str)
    parser.add_argument('--dsnp', help='path of the density SNP file', type=str)
    parser.add_argument('-s', '--start', help='start column number of the founder lines \
                        [defult: --start=4]', default=4, type=int)
    parser.add_argument('-e', '--end', help='end column number of the founder lines \
                        [defult: --end=31]', default=31, type=int)
    parser.add_argument('--header', help='print header or not [default: --header=no]', type=str,
                        choices=['yes','no'], default='no')

    parser.add_argument('-o', '--output', help='output files, like chr1_merged', type=str)
    parser.add_argument('-m','--mode',
                        help='''[default --mode=0];
                        0, for GenSel;
                        1, for PLINK;
                        2, for SNPTEST;
                        3, for raw;
                        ''',
                        choices=[0,1,2,3], default=0, type=int)
    return parser
    #parser = get_parser()
    #parser.print_help()


##########################################################################################
if __name__ == '__main__':
    parser = get_parser()
    args = vars(parser.parse_args())

    if args['diallel'] or args['dsnp'] is not None:
        print(version())
    if args['path'] is not None:
        os.chdir(args['path'])

    ##### read in the pedigree file ######
    st = timeit.default_timer()
    prob=[]
    print("Reading diallel file ...")
    ped = read_ped(pedfile=args['diallel'])
    print("[ ", len(ped), " ] F1 of diallels loaded.")
    write_prob_snp()
    print("[ ", len(prob), " ] SNPs with no polymorphism: write to [ ", outprob,"]")

    ##### read in the density snp file ######
    start=args['start']
    end=args['end']
    print("Reading density SNP file ...")
    dsnp = read2lofd(file_name=args['dsnp'])
    print("[ ", len(dsnp), " ] SNPs loaded from [ ", end-start, " ] founders!")

    ##### start imputation ######
    print("Start SNP imputation")
    output = args['output']
    impute_write(ped=ped, dsnp=dsnp, mode=args['mode'])
    print("File output to: [ ", output, " ]")

    ### test
    #output = "plink"
    #impute_write(ped=ped, dsnp=dsnp[0:10000], mode=1)

    ##### output corresponding phenotype file ########
    write_pheno(mode=args['mode'])
    print("Corresponding phenotype file: [ ", outpheno, " ]")

    #test
    #output ="test"
    #write_pheno(mode=2)

    et = timeit.default_timer()

    print("[ ", "%.0f" % (et - st)/60, " ] minutes of run time!")
    print("imputation finsihed!")



