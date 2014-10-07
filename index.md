---
layout: default
title: Plant3G
subtitle: plant breeding toolset
---

## Overview

The **Plant3G** package was designed to be a transparent engine for SNP manipulation, NGS-enabled genomic selection and genome-wide association studies.

## Motivation



## Features

The ideas are borrowed from other packages, and some of them are
re-implemented in a different way. A selected list of features
include:


### Basic usage and data formats
   - `DSF`: density SNP format (self defined)
   - `VCF`: variant call format
   - `Others`: HapMap format or any other format ...

 
### SNP imputation
   - `Diallel`: from parental high density SNP panel to F1 hybrids. support several downstream analysis software, such as PLINK, SNPTEST, GenSel, ... 
   - `RIL`: from parental high density SNP panel to RILs based on tagging SNPs of RILs.
   - `F1BCn`: from parental high density SNP panel to progenies of the F1BCn based on GBS panel of F1BCn individuals (not implemented yet)

    
### SNP manipulation
   - `MAF`: calculate minor allele frequency
   - `merging`: merging SNPs from different data sets
   - `missing rate`: calculate SNP missing rate for each locus
   - `LD`: calculate linkage disequilibrium (not implemented yet)
   - `Fst`: calculate Fst (not implemented yet)
   - `Others`: ...


### GWAS and GS
3rd party software will be to be install for conducting GWAS and GS. This package will generate a ready-to-run pipeline with user defined options.


   - `single-marker`: single-marker GWAS (**SNPTEST** reuired)
   - `bayes`: bayesian-based GWAS, including bayesA, bayesB and bayesC, each program with many options (**GenSel** required)
   - `GS-training`: training GS model (**GenSel** required)
   - `GS-prediction`: predicting phenotypic performance (**GenSel** required)
   

### Visualization
A collection of visualization tools coded with R.


   - `Manhattan plot`: manhattan plot of GWAS results
   - `locuszoom`: zoom into certain region
   - `Others`: such as regression plot, scatter plot, ...

