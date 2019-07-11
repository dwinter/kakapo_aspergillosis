# Kākāpō _Aspergillus_ analyses

## Installation

Most of the software used in these analyses can be installed with anaconda. 

```sh
conda install -n aspergillus --file conda_pkg.list
```

There are two execptions, picard (the right version is not available in
bioconda) and gatk (you must be registered to download a copy of the .jar).

The Makefile in this directory will install picard. Running `make all`
will handle this, or you can run just this command:

```sh
make include/picard.jar
```

To download `gatk` you will have to register and [go to this
site](https://software.broadinstitute.org/gatk/download/archive) to find version
3.6. You can then run `gatk-register` to provide the location of the .jar file.

## Producing SNPs and assemblies

The SNP calling pipeline replicates the one used by Abdolrasouli A et al. 2015. 
Genomic Context of Azole Resistance Mutations in _Aspergillus fumigatus_ 
Determined Using Whole-Genome Sequencing _mBio_. **6**. doi: [10.1128/mbio.00536-15](dx.doi.org/10.1128/mBio.00536-15).

In short, bwa is used to align reads to a reference genome (part of this repo),
and gatk HaplotypeCaller is used to produce SNPs. The entire analysis can be run
with a single command (this will used 10 CPUs and take a considerable amount of
time)


```sh
make all
```

Steps can be run individually with the Makefile, or with the scripts present in
the `scripts/` directory. 

