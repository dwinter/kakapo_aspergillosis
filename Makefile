REF := ref/Aspergillus_fumigatus.CADRE.12.dna.toplevel.fa
NPROC := 10
# set up names for all reference indices
DICT = $(REF:.fa=.dict)
BWA_IDX = $(addsuffix .bwt, $(REF))
ST_IDX = $(addsuffix .fai, $(REF))

# set up lists of samples, .bams and final and de novo assemblies
SAMPLE_LIST := $(shell cat sample_info)
BAM=$(addsuffix .bam, $(addprefix bam/, $(SAMPLE_LIST)))
VARS=$(addsuffix _filtered_snps_final.vcf.gz, $(addprefix vars/, $(SAMPLE_LIST)))
VARS_GZ=$(addsuffix $(VARS), .gz)
ASSEMBLIES=$(addsuffix _spades.fna, $(addprefix de_novo/, $(SAMPLE_LIST)))


#PREP for analyses
# fetch right picard, set up reference dictionaries
include/picard.jar:
	mkdir -p include
	wget --directory-prefix=include/  https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar

$(DICT): include/picard.jar
	java -jar include/picard.jar CreateSequenceDictionary R=$(REF) O=$(DICT)

$(BWA_IDX):
	bwa index $(REF)

$(ST_IDX):
	samtools faidx $(REF)

.PHONY: indices

indices: $(REF)
	$(MAKE) $(DICT)
	$(MAKE) $(BWA_IDX)
	$(MAKE) $(ST_IDX)

## Perform first alignment 
#
# Not sure if there is a better way to handle multi-target recipies. Workaround
# here is to create a recipe (bams.list) that depends on all bam files
#
bam/%.bam: fq/%_1.fastq.gz fq/%_2.fastq.gz include/picard.jar $(DICT) $(BWA_IDX) $(ST_IDX)
	scripts/make_bam.sh $* $(REF) $(NPROC) > logs/$(*)_ali.log 

bams.list: $(BAM)
	find bam/ -maxdepth 1 -iname *.bam > bams.list
	
.PHONY: alignments

alignments:
	$(MAKE) bams.list



## Variant calling
#
# This recpe produces many files, recipe is written just for the 'final' SNPs

vars/%_filtered_snps_final.vcf.gz: bam/%.bam
	scripts/call_vars.sh bam/$(*).bam $(REF) $(NPROC) > logs/$(*)_vars.log

SNPS.list: $(VARS)
	find vars/ -maxdepth 1 -iname *snps_final.vcf > SNPS.list

.PHONY: vars

vars:
	$(MAKE) SNPS.list

## de novo assembly

de_novo/%_spades.fna: fq/%_1.fastq.gz fq/%_2.fastq.gz
	scripts/assemble.sh $* $(REF) $(NPROC) > logs/$(*)_spades.log

assemblies.list: $(ASSEMBLIES)
	find de_novo/ -maxdepth 1 -iname *_spades.fna > assemblies.list

.PHONY: assemble

assemble:
	$(MAKE) assemblies.list

.PHONY: all

all:
	$(MAKE) alignments
	$(MAKE) vars
	$(MAKE) assemble
	$(MAKE) vars/merged_SNPS.vcf

vars/merged_SNPS.vcf: $(VARS)
	vcf-merge -R "0" $(VARS) > vars/merged_SNPS.vcf


popgen/STRUCTURE.done:
	$(MAKE) vars/merged_SNPS.vcf
	# This gets a bit hacky...
	#really need to work out why pyVCF is failing in python3...
	python2 scripts/vcf_to_STRUCTURE.py popgen/merged
	cd popgen/
	seq 2 10 | parallel -j $(NRPOC) ../scripts/run_STRUCTURE_pa.sh {}


clean_tempfiles:
	rm  vars/temp/*
	rm  bam/temp/*
