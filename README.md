# Kākāpō _Aspergillus_ analyses

## Installation

Most of the software used in these analyses can be installed with anaconda. 

```sh
conda install -n aspergillus --file conda_pkg.list
```

There are two execptions, picard (the right version is not available in
bioconda) and gatk (you must be registered to download a copy of the .jar).

The Makefile in this directory will install picard (indeed running `make all`
will handle this).

```sh
make include/picard.jar
```

To download gatk you will have to register and [go to this
site](https://software.broadinstitute.org/gatk/download/archive) to find version
3.6. You can then run `gatk-register` to provide the location of the .jar file.

## Producing SNPs and assemblies

```sh
make all
```
