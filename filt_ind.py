#!/usr/bin/env python

'''
Extract the genotypes of a list of individuals from a Dosafe files having
IDs as colnames and Markers in each row
Roberto Lozano ... rjl278@cornell.edu .. for bugs
'''

import re
import sys
from operator import itemgetter

# Arguments:
lista  = sys.argv[1]      # Your list with the desired individuals

quieroesto = list()
with open(sys.argv[1]) as f:
    for lines in f:
        a = lines.split("\t")
        quieroesto.append(a[0].rstrip())

indices = list()
with open("/home/NVME/Pro-seq/Quantitative/Predictions/genos/IDs.imputed") as g:
    #Set a counter:
    n = 1    
    for lines in g:
        a = lines.split("\t")
        if a[0].rstrip() in quieroesto :
            indices.append(n)
            #print (n)
        n += 1

indices.insert(0, 0)

with open("/home/NVME/Pro-seq/Quantitative/Predictions/genos/C1.I2.dosage") as f:    
    for lines in f:
        imprime = list()
        a = lines.split("\t")
        n = 0
        for elements in a:
            if n in indices:
                imprime.append(elements.rstrip())    
            n +=1
        print ("\t".join(imprime))      
