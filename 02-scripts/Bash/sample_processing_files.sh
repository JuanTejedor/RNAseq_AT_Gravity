#!/bin/bash
# Estos parametros hay que añadirlos al lanzar el scripts en el sevidor

# --job-name=s1
# --job-output=s1
# se lanza y los parámetros que requiera el script

#SBATCH --export=ALL

## This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International$
## Author: Francisco J. Romero-Campero
## Date: May 2023
## Email: fran@us.es

SAMPLE=$1
FICHERO_R1=$2
FICHERO_R2=$3

echo "Sample" $SAMPLE
echo /home/omicas/grupo6/samples/sample$SAMPLE

echo "fichero1"
echo $FICHERO_R1

echo "fichero2"
echo $FICHERO_R2


## Downloading sample file
cd /home/omicas/grupo6/samples/sample$SAMPLE

echo "He accedido a la carpeta"

## Para descargar en SRA  # En nuestro paper es con enlace wget
# De todas formas lo dejo apuntado
# fastq-dump --gzip --split-files <accession_number_SRA el que pone RUN>

## Para descargar en otra base de datos # Este es nuestro caso
mv $FICHERO_R1 R1.fastq.gz

echo "file 1 pasado a R1.fastq.gz"

mv $FICHERO_R2 R2.fastq.gz

echo "file 2 pasado a R1.fastq.gz"

## Sample quality control and read mapping to reference genome
if [ -f R1.fastq.gz ]
then
   fastqc R1.fastq.gz
   fastqc R2.fastq.gz
# Hay que modificar esta parte para nuestro script
   hisat2 --dta -x ../../genome/index -1 R1.fastq.gz -2 R2.fastq.gz -S sample$SAMPLE.sam
fi

## Generting sorted bam file
samtools sort --threads 5 -o sample$SAMPLE.bam sample$SAMPLE.sam
rm sample$SAMPLE.sam
rm *.fastq.gz
samtools index sample$SAMPLE.bam
#bamCoverage -bs 10 --normalizeUsing CPM --bam sample<N>.bam -o sample<N>.bw


## Transcript assembly
stringtie -G ../../annotation/annotation.gtf -o sample$SAMPLE.gtf -l sample$SAMPLE sample$SAMPLE.bam

## Gene Expression Quantification
stringtie -e -B -G ../../annotation/annotation.gtf -o sample$SAMPLE.gtf sample$SAMPLE.bam

rm *.bam
