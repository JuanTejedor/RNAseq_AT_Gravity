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


SAMPLE=$1 # NUMERO ID DE LA MUESTRA
ENLACES=$2 # FICHERO QUE CONTIENE LOS 2 ENLACES: R1 Y R2



SAMPLE_DIR=$(find ~/ -type d -name 'samples' | grep 'grupo6' | grep 'rna')/sample$SAMPLE



echo "Muestra: $SAMPLE"
echo $SAMPLE_DIR


cd $SAMPLE_DIR
echo "El directorio actual es $pwd"


## Downloading sample file
wget --tries=10 -i $ENLACES

FICHERO_R1=$(ls *R1*) # NOMBRE DE FICHERO R1
FICHERO_R2=$(ls *R2*) # NOMBRE DEL FICHERO R2

echo "Fichero R1: $FICHERO_R1"
echo "Fichero R2: $FICHERO_R2"
## Para descargar en SRA  # En nuestro paper es con enlace wget
# De todas formas lo dejo apuntado
# fastq-dump --gzip --split-files <accession_number_SRA el que pone RUN>

## Para descargar en otra base de datos # Este es nuestro caso
mv $FICHERO_R1 R1.fastq.gz

echo "$FICHERO_R1 renombrado como R1.fastq.gz"

mv $FICHERO_R2 R2.fastq.gz

echo "$FICHERO_R2 renombrado como R2.fastq.gz"

## Sample quality control and read mapping to reference genome
if [ -f R1.fastq.gz ]
then
	echo "------------------------------------------------------------------------------------------------------------------------------------------"
	echo "##################### Se ejecuta FASTQC #####################"
	echo "------------------------------------------------------------------------------------------------------------------------------------------"
	fastqc R1.fastq.gz
	fastqc R2.fastq.gz
# Hay que modificar esta parte para nuestro script
	echo "------------------------------------------------------------------------------------------------------------------------------------------"
	echo "##################### Se ejecuta HISAT2 #####################"
	echo "------------------------------------------------------------------------------------------------------------------------------------------"
	hisat2 --dta -x ../../genome/index -1 R1.fastq.gz -2 R2.fastq.gz -S sample$SAMPLE.sam -p 6
fi

## Generting sorted bam file
echo "------------------------------------------------------------------------------------------------------------------------------------------"
echo "##################### Se ejecuta SAMTOOLS (sam a bam) #####################"
echo "------------------------------------------------------------------------------------------------------------------------------------------"
samtools sort --threads 6 -o sample$SAMPLE.bam sample$SAMPLE.sam
echo "------------------------------------------------------------------------------------------------------------------------------------------"
echo "##################### Se eliminan .sam y .fastq.gz #####################"
echo "------------------------------------------------------------------------------------------------------------------------------------------"
rm sample$SAMPLE.sam
rm *.fastq.gz
echo "------------------------------------------------------------------------------------------------------------------------------------------"
echo "##################### Se ejecuta SAMTOOLS (indice) #####################"
echo "------------------------------------------------------------------------------------------------------------------------------------------"
samtools index sample$SAMPLE.bam
#bamCoverage -bs 10 --normalizeUsing CPM --bam sample<N>.bam -o sample<N>.bw

## Transcript assembly
echo "------------------------------------------------------------------------------------------------------------------------------------------"
echo "##################### Se ejecuta STRINGTIE (ensamblaje) #####################"
echo "------------------------------------------------------------------------------------------------------------------------------------------"
stringtie -G ../../annotation/annotation.gtf -o sample$SAMPLE.gtf -l sample$SAMPLE sample$SAMPLE.bam

## Gene Expression Quantification
echo "------------------------------------------------------------------------------------------------------------------------------------------"
echo "##################### Se ejecuta STRINGTIE (cuantificacion) #####################"
echo "------------------------------------------------------------------------------------------------------------------------------------------"
stringtie -e -B -G ../../annotation/annotation.gtf -o sample$SAMPLE.gtf sample$SAMPLE.bam

echo "------------------------------------------------------------------------------------------------------------------------------------------"
echo "##################### Se eliminan .bam #####################"
echo "------------------------------------------------------------------------------------------------------------------------------------------"
rm *.bam
