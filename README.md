Custom scripts for the NGS analysis, pausing index calculations, metagene plots and genomic partitioning of the working article titled:

# RNA polymerase mapping in plants identifies enhancers enriched in causal variants

## NGS Pipeline:
Pipeline_cassava.sh

## Metagene Plots:


## Methylation analysis:


## Genomic Partitioning:

random.py
Given a bed file with some regions will generate relationship matrices using LDAK as described in the figure S8

Genomic_partitioning.R
Will calculate Vu, Ve, and the weights for the environmental effects, focal kernel and rest of the genome kernel.
use as: Rscript Genomic_partitioning.R <focal grm> <rest of the genome grm> <TRAIT> <output>
  
phenos.Rdata
Phenotypes present on the supplementary table S2
