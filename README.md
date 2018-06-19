Custom scripts for the NGS analysis, pausing index calculations, metagene plots and genomic partitioning of the working article titled:

# RNA polymerase mapping in plants identifies enhancers enriched in causal variants

## NGS Pipeline:
**Pipeline_cassava.sh:**
Script used to proccess and align the raw sequencing libraries. See *Analysis of NGS data* for details.

## Main text Figures:

**Figure1.Rmd:** R notebook code used to generate Fig 1 of the manuscript

**Figure2.Rmd:** R notebook code used to generate Fig 2 of the manuscript

**Figure3.Rmd:** R notebook code used to generate Fig 3 of the manuscript

## Supplementary Figures:

**Supplementary_figures.Rmd:** R notebook code used to generate supplementary figures 1 - 6

## Genomic Partitioning:

**random.py:** 
Given a bed file with some regions will generate relationship matrices using LDAK as described in the figure S8

**Genomic_partitioning.R:** 
Will calculate Vu, Ve, and the weights for the environmental effects, focal kernel and rest of the genome kernel.
use as: `Rscript Genomic_partitioning.R <focal grm> <rest of the genome grm> <TRAIT> <output>`
  
**phenos.Rdata:** 
Phenotypes present on the supplementary table S2

## Misc:
**GERP_formatting.py:**
Format the original GERP files to location aware bed files, only showing position with GERP scores > 1

**Pausing_index.py:**
This Code will produce Pausing index and divergent index scores. 

**Genomic_Partitioning_protocol.Rmd:**
Old version of the Genomic Partition pipeline

**filt_ind.py:**
Extract the genotypes of a list of individuals from a Dosage file having IDs as "colnames" and Markers as "rownames"

