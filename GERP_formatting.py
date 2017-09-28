#!/usr/bin/env python

'''
Format the original GERP files to location aware bed files
Only showing position with GERP scores > 1
Original files can be find @ Tuber: "/export/species/Manihot_esculenta_old/gbs/IGDbuildNewWithV6/GERPscore_cassavaV6"
Roberto Lozano ... rjl278@cornell.edu
'''

#Load some packages
import re
import sys

#Arguments
gerp       = sys.argv[1]    # GERP file
chromosome = sys.argv[2]    # Chromosome being evaluated 


n = 1 # counter
with open(gerp) as h:
    for lines in h:
        a = lines.split("\t")
        if float(a[1]) > 0 :
            name = "S"+chromosome+"_"+str(n)
            print(chromosome, n, n+1, name, a[1].rstrip(), "+", sep = "\t", end = "\n")
        n +=1