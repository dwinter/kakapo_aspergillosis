#!/bin/sh
#
# Perform a quick assembly for a set of paired end reads with spades
# Usage:
#   assemble.sh [sample-name] [reference-genome] [nproc]
# Where sample name is the stem of the .fastq file in /fq
#       reference-genome is a high-quality reference used for comparisons
#       nproc = number of processors to use
#

SAMPLE=$1
REF=$2
NPROC=$3
R1=fq/${SAMPLE}_1.fastq
R2=fq/${SAMPLE}_2.fastq

spades.py -t 10 -m 24 -k  33,55,77,99,127 --careful --pe1-1 ${R1} --pe1-2 ${R2} -o de_novo/${SAMPLE}
quast -t 10 -r ${REF} -o de_novo/${SAMPLE}/quast  de_novo/${SAMPLE}/scaffolds.fasta
scripts/discard_chaff.py de_novo/${SAMPLE}/scaffolds.fasta 500 > de_novo/${SAMPLE}_spades.fna
