#!/usr/bin/env python

'''
This Code will produce Pausing index scores
The indexes will be calculated for both strands
Adding new changes
Checking rsa key

'''


import re
import os
import sys
from operator import itemgetter
import subprocess
import HTSeq

# Arguments:
lista  = sys.argv[1]      # Bed file with each gene 

## sudo ln -fs /usr/lib/libcurl.so.4 /usr/local/lib    => This is to prevent error messages cause I have two curl installs

##############################################################
#I. CALCULATE THE NUMBER OF EXONS AND CDS LENGTH USING HTSEQ
##############################################################

gff_file = HTSeq.GFF_Reader("/home/DB2/RNAseq_V2/DB/mesculenta_305_v6.1.gene_exons.gff3")

transcripts = {}

for feature in gff_file:
   if feature.type == "exon":
      transcript_id = feature.attr['Parent']
      if transcript_id not in transcripts:
         transcripts[ transcript_id ] = list()
      transcripts[ transcript_id ].append( feature )
      
distancias = dict()

for transcript_id in sorted( transcripts ):      
   transcript_length = 0
   for exon in transcripts[ transcript_id ]:
      transcript_length += exon.iv.length + 1
   distancias[transcript_id] =  transcript_length , len( transcripts[ transcript_id ] )
   #print (transcript_id, transcript_length , len( transcripts[ transcript_id ]))

#Removing the alternative transcripts and shortening the gene names to standardize
#The dictionary containing the length of the CDS and the exon Length = dicready

dicready = dict()
for keys in distancias.keys():
   if re.search(".1.v6.1", keys):
      dicready[keys[0:15]] = distancias[keys][0] , distancias[keys][1]


##############################################################
#II. CALCULATE PAUSING INDEX AND DIVERGENT INDEX
##############################################################

header = "Gene", "Pause", "Body", "PI", "Divergent", "NonDivergent", "DTI", "Expression", "GeneLength", "mRNALength", "nexons",  "Strand"
print("\t".join(header))

with open(sys.argv[1]) as f:
    for lines in f:
        a = lines.split("\t")
        length = int(a[2]) - int(a[1])
        
        # Calculate reads mapping to each gene and divide that by the lenght of the gene
        Coor     = a[0]+ ":" + str(int(a[1])) + "-" + str(int(a[2]))
        nreads = os.popen("samtools6 view -c white_trimmed_sorted.bam %s" % (Coor)).read().rstrip()
        expression =  "%.2f" % round(int(nreads)/(int(length)/1000), 2)
        
        
        ## POSITIVE STRAND
        
        if a[4].rstrip() == "+" and int(a[1])>800 and int(length) > 300:
            
            #CALCULATING PAUSING INDEX
            pausing     = a[0]+ ":" + str(int(a[1])-100) + "-" + str(int(a[1]) + 300)
            pauslen     = (int(a[1])+300) - (int(a[1])-100)
            
            body        = a[0]+ ":" + str(int(a[1])+300) + "-" + str(int(a[2]))    
            bodylen     = int(a[2]) - (int(a[1])+300)
            
            #Ccoverages
            pausingd = os.popen("samtools6 depth -a -r %s white_plus.bam | awk '{sum += $3;print sum}' |tail -n 1" % (pausing)).read().rstrip()
            if pausingd.rstrip() == "":
                pausingd = 0
            bodyd     = os.popen("samtools6 depth -a -r %s white_plus.bam | awk '{sum += $3;print sum}' |tail -n 1" % (body)).read().rstrip()
            if bodyd == "":
                bodyd = 0
            
            try:
                index =    str(round(     (int(pausingd)/pauslen)/(int(bodyd)/bodylen), 3))
            except ZeroDivisionError:
                index = -1
            
            
            #CALCULATING DIVERGENT INDEX
            Ndivergent       = a[0]+ ":" + str(int(a[1])-300) + "-" + str(int(a[1]) + 300)
            divergent        = a[0]+ ":" + str(int(a[1])-1000) + "-" + str(int(a[1])) 
            
            #calculate the coverage for each segment
            Ndivergentd = os.popen("samtools6 depth -a -r %s white_plus.bam | awk '{sum += $3;print sum}' |tail -n 1" % (Ndivergent)).read().rstrip()
            if Ndivergentd.rstrip() == "":
                Ndivergentd = 0
            
            divergentd     = os.popen("samtools6 depth -a -r %s white_minus.bam | awk '{sum += $3;print sum}' |tail -n 1" % (divergent)).read().rstrip()
            if divergentd == "":
                divergentd = 0
            
            try:
                indexD =    str(round(     (int(divergentd))/(int(Ndivergentd)), 3))
            except ZeroDivisionError:
                indexD = -1
            
            
            newline = (a[3], str(pausingd), str(bodyd), str(index),  str(divergentd), str(Ndivergentd), str(indexD), str(expression), str(length), str(dicready[a[3]][0]), str(dicready[a[3]][1]),  a[4].rstrip() )
            print("\t".join(newline))
        
        
        ## NEGATIVE STRAND
            
        if a[4].rstrip() == "-" and int(a[1])>800 and int(length) > 300:
            #Calculating Pausing coordinates for +
            pausing     = a[0]+ ":" + str(int(a[2])-300) + "-" + str(int(a[2]) + 100)
            pauslen     = (int(a[2])+100) - (int(a[2])-300)
           
            body        = a[0]+ ":" + str(int(a[1])) + "-" + str(int(a[2])-300)    
            bodylen     = (int(a[2])-300) - int(a[1])
            
            #calculate the coverage for each segment
            pausingd = os.popen("samtools6 depth -a -r %s white_minus.bam | awk '{sum += $3;print sum}' |tail -n 1" % (pausing)).read().rstrip()
            if pausingd.rstrip() == "":
                pausingd = 0
            bodyd     = os.popen("samtools6 depth -a -r %s white_minus.bam | awk '{sum += $3;print sum}' |tail -n 1" % (body)).read().rstrip()
            if bodyd == "":
                bodyd = 0
            
            try:
                index =    str(round(     (int(pausingd)/pauslen)/(int(bodyd)/bodylen), 3))
            except ZeroDivisionError:
                index = -1
            
            #CALCULATING DIVERGENT INDEX
            Ndivergent       = a[0]+ ":" + str(int(a[2])-300) + "-" + str(int(a[2]) + 300)
            divergent        = a[0]+ ":" + str(int(a[2])) + "-" + str(int(a[2])+1000) 
            
            #calculate the coverage for each segment
            Ndivergentd = os.popen("samtools6 depth -a -r %s white_minus.bam | awk '{sum += $3;print sum}' |tail -n 1" % (Ndivergent)).read().rstrip()
            if Ndivergentd.rstrip() == "":
                Ndivergentd = 0
            
            divergentd     = os.popen("samtools6 depth -a -r %s white_plus.bam | awk '{sum += $3;print sum}' |tail -n 1" % (divergent)).read().rstrip()
            if divergentd == "":
                divergentd = 0
            
            try:
                indexD =    str(round(     (int(divergentd))/(int(Ndivergentd)), 3))
            except ZeroDivisionError:
                indexD = -1
            
            
            newline = (a[3], str(pausingd), str(bodyd), str(index),  str(divergentd), str(Ndivergentd), str(indexD), str(expression), str(length), str(dicready[a[3]][0]), str(dicready[a[3]][1]), a[4].rstrip() )
            print("\t".join(newline))

        
        
        
        
        
        
        
        