---
layout: default
title: Objects
subtitle: SNP and InDel Objects to manipulate
---

For most of the plant genomic and genetic study, a simple SNP matrix with SNP ID and SNP variations for multiple individuals is sufficient enough. For this purpose, **Plant3G** defined a simple format (**BED+**) for SNP manipulation. This format borrowed the simplicity and flexibility of [BED format](http://bedtools.readthedocs.org/en/latest/content/general-usage.html). We provide some tools to transform from and into other formats.

Below is the description of the **BED+ format**. Modified from [UCSC](http://genome.ucsc.edu/FAQ/FAQformat#format1).

1. **chr** - The name of the chromosome (e.g. chr1, 1). 
2. **start** - The starting position (zero-based) of the feature in the chromosome.
- The first base in a chromosome is numbered 0.
- For example, start=9, end=20 is interpreted to span bases 10 through 20, inclusive. 
3. **end** - The one-based ending position of the feature in the chromosome.
4. **ID** - Defines the name of the BED feature.
5. **+** - Unlimited columns of plants with Genotype information. Each column represent one plant.

-----














