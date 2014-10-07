zmSNPtools
==========

A collection of python codes designed to process a range of basic, \
large-scale analyses for SNPs of maize diversity panel.

Usage:
=========
Copy the piece of python code to the directory in your $PATH:
for example:
```
$: cp impute4diallel/impute4diallel.py ~/bin/impute4diallel

```
Then, you will be ready to run it
```
$: impute4diallel -help
```

Description
=========
### fpSNP
find a set of population specific fingerprint SNPs

## update of the program
10/3/2014: v0.1 python code tested with a small set of simulated data

## usage:

```
git clone git@github.com:yangjl/zmSNPtools.git
cp zmSNPtools/packages/fpSNP/fpSNP_v0.1.py ~/bin/fpSNP
export PATH=$PATH:~/bin/
fpSNP -h

```



### impute4diallel
SNP Imputation for maize diallel population from parental SNP panel to F1 hybrids.

### dsf2GWAS
transform density SNP format (DSF) to GWAS (gensel4.2 supported) and other formats

### 




Update
=========
### General
1. changed gitignore
2. unwatch the MyPack and test folders


### packages:
1. updated a SNP merge package **snp3merge**: merge hapmap1, hapmap2 and RNA-seq data
2. updated **impute4diallel**
3. updated **snpfrq**
4. updated **merge4ames**

### MyPack
1. updated and tested **snp_sampling** program for XBSA project;


