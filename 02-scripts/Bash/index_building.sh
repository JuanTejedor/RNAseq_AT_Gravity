#!/bin/bash
#SBATCH --job-name=index_grupo6
#SBATCH --export=ALL
#SBATCH --output=genome_index_grupo6

## This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International$
## Author: Francisco J. Romero-Campero
## Date: June 2021
## Email: fran@us.es

## Dowloading reference genome and annotation
cd /home/omicas/grupo6/genome 
wget -O genome.fa.gz https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip genome.fa.gz

cd /home/omicas/grupo6/annotation
wget -O annotation.gtf.gz https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.58.gtf.gz 
gunzip annotation.gtf.gz

## Building reference genome index
cd /home/omicas/grupo6/genome
extract_splice_sites.py ../annotation/annotation.gtf > splices_sites.ss
extract_exons.py ../annotation/annotation.gtf > exons.exon
hisat2-build --ss splices_sites.ss --exon exons.exon genome.fa index
