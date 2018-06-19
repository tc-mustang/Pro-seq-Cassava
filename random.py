#!/usr/bin/env python

import sys
import os

os.system("ldak5 --cut-genes random --bfile /home/NVME/Pro-seq/Quantitative/dREG/LD/Enhancers/Enhancer --genefile random.ldak --ignore-weights YES")
snp_set = "random/genes.predictors.used" 

## Calculate a set of markers that are in LD with the random set
os.system("plink --bfile /home/NVME/Pro-seq/Quantitative/dREG/LD/GGC1 --show-tags %s --tag-kb 100 --tag-r2 0.9" % snp_set)
os.system("cat %s plink.tags | sort | uniq -u > exclude.snps" % snp_set)
os.system("cat exclude.snps /home/NVME/Pro-seq/Quantitative/dREG/LD/Random/snpsingenes/kernel1.must | sort | uniq -d > tmp && cat tmp exclude.snps |sort |uniq -u > exclude.this")

## Make plink files withouth this markers
os.system("plink --bfile /home/NVME/Pro-seq/Quantitative/dREG/LD/GGC1 --exclude exclude.this --make-bed --out rplink")

##Generate the LDAK weigths!!!
os.system("awk '{print $2}' rplink.bim > SNPs_rplink")
os.system("ldak5 --cut-weights sectionGP --bfile rplink --extract SNPs_rplink")

sections = os.popen("more sectionGP/section.number" ).read().rstrip()

os.system("for j in {1..%s}; do ldak5 --calc-weights sectionGP/ --bfile rplink --extract SNPs_rplink --section $j; done" % sections)
os.system("ldak5 --join-weights sectionGP --bfile rplink --extract SNPs_rplink")

## Calculate relationship Matrices
os.system("ln -s random/genes.predictors.used ./kernel2 ") 
os.system("cat SNPs_rplink kernel2 | sort | uniq -u > kernel1 ")

os.system("ldak5 --cut-kins partitionGP --bfile rplink --partition-number 2 --partition-prefix kernel")
os.system("ldak5 --calc-kins partitionGP --bfile rplink --weights sectionGP/weights.all --partition 1  --kinship-raw YES --power -0.25")
os.system("ldak5 --calc-kins partitionGP --bfile rplink --weights sectionGP/weights.all --partition 2  --kinship-raw YES --power -0.25")

