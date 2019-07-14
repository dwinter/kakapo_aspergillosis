#!/bin/sh 
#
# Align readsfrom sample (first argument) to reference genome (second argumetn).
#
# Usage:
#   make_bam.sh [sample_name] [reference genome] [nproc]
#
# This scrip asssumes paired end reads with names fq/[sample_name]_1.fq and
# fq/[sample_name]_2.fq exist and that the directory with a refernce genome
# includes bwa, samtools and picard indexes of the genome
#
# Creates a file bam/[sample-name].bam


SAMPLE=$1
REF=$2
NPROC=$3
R1=fq/${SAMPLE}_1.fastq.gz
R2=fq/${SAMPLE}_2.fastq.gz
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:Illumina\tLB:${SAMPLE}"
#Align
echo bwa mem -M -t ${NPROC} -R "${RG}" ${REF} ${R1} ${R2} 
bwa mem -M -t ${NPROC} -R "${RG}" ${REF} ${R1} ${R2} | samtools sort -"@"${NPROC} -T /tmp/tmp.aln -o bam/temp/${SAMPLE}_raw.bam
samtools index bam/temp/${SAMPLE}_raw.bam
#clean up, first duplicates
java -jar include/picard.jar MarkDuplicates INPUT=bam/temp/${SAMPLE}_raw.bam OUTPUT=bam/temp/${SAMPLE}_MD.bam METRICS_FILE=stats/${SAMPLE}_MD_info.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=TRUE
#..then realign indels
gatk -T RealignerTargetCreator -nt ${NPROC} -R ${REF} -I bam/temp/${SAMPLE}_MD.bam -o stats/${SAMPLE}.realignment_targets.list
gatk -T IndelRealigner  -R ${REF} -I bam/temp/${SAMPLE}_MD.bam -targetIntervals stats/${SAMPLE}.realignment_targets.list -o bam/${SAMPLE}.bam

