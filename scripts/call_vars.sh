#!/bin/sh
#
# Call variants using the pipeline used for existing A. fumigatus diversity 
# databases. 
#
# Usage:
#   call_vars.sh [bam-file] [reference-genome] [nproc]
#
# the bam file name should have the form "sample-name.bam" and .vcf files
# will be written to "vars/" directory

BAM=$1
REF=$2
NPROC=$3
SAMPLE=`basename ${BAM} .bam`
RECAL_BAM=bam/${SAMPLE}_recal.bam 

echo "=== Initial SNP calling for ${SAMPLE} ====" 1>&2
gatk -T HaplotypeCaller -nct ${NPROC} -R ${REF} -I ${BAM} -o vars/temp/${SAMPLE}_raw_variants.vcf -pcrModel NONE -ploidy 1 -stand_call_conf 30 -mbq 20 -A QualByDepth -A AlleleBalance -XL ref/Aspergillus_fumigatus.CADRE.12.dna.toplevel.repeat.intervals

#Extract SNPs and INDELs
gatk -T SelectVariants -R ${REF} -V vars/temp/${SAMPLE}_raw_variants.vcf -selectType SNP -o vars/temp/${SAMPLE}_raw_snps.vcf
gatk -T SelectVariants -R ${REF} -V vars/temp/${SAMPLE}_raw_variants.vcf -selectType INDEL -o vars/temp/${SAMPLE}_raw_indels.vcf

gatk -T VariantFiltration -R ${REF} -V vars/temp/${SAMPLE}_raw_snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0" --filterName LowConf -o vars/temp/${SAMPLE}_filtered_snps.vcf

gatk -T VariantFiltration -R ${REF} -V vars/temp/${SAMPLE}_raw_indels.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filterName LowConf -o vars/temp/${SAMPLE}_filtered_indels.vcf

echo "=== BSQR for ${SAMPLE} ====" 1>&2
#Base Quality Score Recalibration (BQSR) #1
gatk -T BaseRecalibrator -R ${REF} -I ${BAM} -knownSites vars/temp/${SAMPLE}_filtered_snps.vcf -knownSites vars/temp/${SAMPLE}_filtered_indels.vcf -o vars/temp/${SAMPLE}_recal_data.table

#Base Quality Score Recalibration (BQSR) #2
gatk -T BaseRecalibrator -R ${REF} -I ${BAM} -knownSites vars/temp/${SAMPLE}_filtered_snps.vcf -knownSites vars/temp/${SAMPLE}_filtered_indels.vcf -BQSR vars/temp/${SAMPLE}_recal_data.table -o vars/temp/${SAMPLE}_post_recal_data.table

#Apply BQSR
gatk -T PrintReads -R ${REF} -I ${BAM} -o ${RECAL_BAM} -BQSR vars/temp/${SAMPLE}_recal_data.table


echo "=== Final variant calling for ${SAMPLE} ====" 1>&2
gatk -T HaplotypeCaller -nct ${NPROC} -R ${REF} -I ${RECAL_BAM} -o vars/temp/${SAMPLE}_raw_variants_recal.vcf -pcrModel NONE -ploidy 1 -stand_call_conf 30 -mbq 20 -A QualByDepth -A AlleleBalance -XL ref/Aspergillus_fumigatus.CADRE.12.dna.toplevel.repeat.intervals


gatk -T SelectVariants -R ${REF} -V vars/temp/${SAMPLE}_raw_variants_recal.vcf -selectType SNP -o vars/temp/${SAMPLE}_raw_snps_recal.vcf

gatk -T SelectVariants -R ${REF} -V vars/temp/${SAMPLE}_raw_variants_recal.vcf -selectType INDEL -o vars/temp/${SAMPLE}_raw_indels_recal.vcf

#Annotate SNPs only
gatk -T VariantAnnotator -R ${REF} --variant vars/temp/${SAMPLE}_raw_snps_recal.vcf -o vars/temp/${SAMPLE}_annotated_snps.vcf --annotation AlleleBalance

#Filter SNPs
gatk -T VariantFiltration -R ${REF} -V vars/temp/${SAMPLE}_annotated_snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0 || DP < 10 || ABHom < 0.8" --filterName LowConf -o vars/${SAMPLE}_filtered_snps_final.vcf

#Filter INDELs
gatk -T VariantFiltration -R ${REF} -V vars/temp/${SAMPLE}_raw_indels_recal.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filterName LowConf -o vars/${SAMPLE}_filtered_indels_final.vcf

grep PASS vars/${SAMPLE}_filtered_snps_final.vcf | awk '$4=="A"||$4=="C"||$4=="G"||$4=="T"' | awk '$5=="A"||$5=="C"||$5=="G"||$5=="T"' > vars/${SAMPLE}_final_snps.body
grep "#" vars/${SAMPLE}_filtered_snps_final.vcf > vars/${SAMPLE}_final.head
cat vars/temp/${SAMPLE}_final.head vars/${SAMPLE}_final_snps.body > vars/${SAMPLE}_final_snps.vcf

grep PASS vars/${SAMPLE}_filtered_indels_final.vcf | awk
'$4=="A"||$4=="C"||$4=="G"||$4=="T"' | awk '$5=="A"||$5=="C"||$5=="G"||$5=="T"'> vars/${SAMPLE}_final_indels.body
grep "#" vars/${SAMPLE}_filtered_indels_final.vcf > vars/${SAMPLE}_final_indels.head
cat vars/${SAMPLE}_final_indels.head vars/${SAMPLE}_final_indels.body > vars/${SAMPLE}_final_indels.vcf

#zip and index SNPs
bgzip vars/${SAMPLE}_filtered_snps_final.vcf 
tabix -p vcf varsvars/${SAMPLE}_filtered_snps_final.vcf.gz

gatk -T DepthOfCoverage -I ${RECAL_BAM} -o stats/${sample}_coverage -R ${REF}
samtools flagstat ${RECAL_BAM} > stats/${sample}_flagstat

