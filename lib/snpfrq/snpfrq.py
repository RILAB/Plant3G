#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
import sys
import argparse
import textwrap
import os

def version():
    v1 = """
    ##########################################################################################
    snpfrq version 0.2.0
    Jinliang Yang
    updated: Oct.24.2013
    --------------------------------

    compute SNP frq and missingness
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
    parser.add_argument('-m','--mode',
                        help='''[default --mode=2];
                        0, for GenSel [Not currently support];
                        1, for PLINK;
                        2, for SNPTEST;
                        ''',
                        choices=[1,2], default=2, type=int)
    
    parser.add_argument('-o', '--output', help='output files, like chr1_merged', type=str)
    
    return parser
    #parser = get_parser()
    #parser.print_help()

def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

##########################################################################################

def read2lofd(file_name):
    infile = open(file_name, 'r')
    
    firstline = infile.readline()
    first = firstline.split()
    #check the header information#
    secline = infile.readline()
    sec = secline.split()
    if(len(str(sec[4])) != 1):
        warning("SNP should code with one latter, like:'A'!")
    infile.close()
    
    temp = []
    infile = open(file_name, 'r')
    next(infile) # skip the header line
    for line in infile:
        tokens = line.split()
        out = get_loci_info(tokens)
        if out:
            temp.append(out)
    infile.close()
    return temp
#test: test = read2lofd(f)



def read_snptest(file_name):
    infile = open(file_name, 'r')
    fl = infile.readline()
    line = fl.split()
    if (len(line) -5) %3 != 0:
        warning("SNPTEST formatting error!")
    else:
        print("Loading SNP with SNPTEST format!")
    infile.close()
    
    infile=open(file_name, 'r')
    for line in infile:
        line = line[:-1]
        tokens = line.split('\t')
        



def compute_snptest(loci):
    
    miss = nomiss = a1 = a2 =0
    for snp in xrange(5, len(loci)):
        loci[snp] = loci[snp].split()
        if set(loci[snp]) == {0}:
            miss += 1
        else:
            nomiss +=1
            a1 = a1 + 2* float(loci[snp][0]) + float(loci[snp][1])
            a2 = a2 + float(loci[snp][1]) + 2*float(loci[snp][2])
    
    res = {loci[3]:a1, loci[4]:a2, 'missing':miss, 'nomiss':nomiss, 'total':len(loci)-5}
    
    return res

tokens = line.split('\t')    
compute_frq_missing(tokens)    


def compute_plink(loci):
    
    miss = nomiss = a1 = a2 =0
    for snp in xrange(5, len(loci)):
        loci[snp] = loci[snp].split()
        if set(loci[snp]) == {0}:
            miss += 1
        else:
            nomiss +=1
            a1 = a1 + 2* float(loci[snp][0]) + float(loci[snp][1])
            a2 = a2 + float(loci[snp][1]) + 2*float(loci[snp][2])
    
    res = {loci[3]:a1, loci[4]:a2, 'missing':miss, 'nomiss':nomiss, 'total':len(loci)-5}
    
    return res

tokens = line.split('\t')    
compute_frq_missing(tokens)    





def get_loci_info(tokens):
    
    set0 = set(tokens[start:end])
    if 'N' in set0:
        set0.remove('N')
    snpset = list(set0)
        
    temp1 = {} 
    if len(snpset) != 2:
        prob.append(tokens)
    else: 
        c1 = tokens[start:end].count(snpset[0])
        c2 = tokens[start:end].count(snpset[1])
   
        for i in range(0, len(first)):
            temp1[first[i]] = tokens[i]
            if c1 >= c2:
                temp1['major'] = snpset[0]
                temp1['minor'] = snpset[1]
                temp1['maf'] = c2/(c1+c2)
            else:
                temp1['major'] = snpset[1]
                temp1['minor'] = snpset[0]
                temp1['maf'] = c1/(c1+c2)
            temp1['missing'] = end - start - c1 - c2
    return temp1

##########################################################################################

    



    
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




