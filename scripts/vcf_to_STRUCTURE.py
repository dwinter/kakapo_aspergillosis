#!/usr/bin python
import sys
import vcf

"""
Convert a vcf file to STRUCTURE input format. 

NOTE: This script is written to produce data from haploid individuals and has
only been tested on vcf files form GATK HaplotypeCaller. IT IS VERY LIKELY THIS
SCRIPT WILL BREAK OR PRODUCE UNEXPECTED RESULTS WHEN USED WITH OTHER DATA

Usage
vcf_to_STRUCTURE.py [in.vcf] [out-file-stem]

Output:
    Two files one "stem.structure" which contains the unput data, the other 
    "stem.mparams" contains a 'mianparams' paramter file for STRUCTURE

"""

#Default params file (but maxpops = 4, numinds and numloci are templated)

params_template = """

KEY PARAMETERS FOR THE PROGRAM structure.  YOU WILL NEED TO SET THESE
IN ORDER TO RUN THE PROGRAM.  VARIOUS OPTIONS CAN BE ADJUSTED IN THE
FILE extraparams.


"(int)" means that this takes an integer value.
"(B)"   means that this variable is Boolean 
        (ie insert 1 for True, and 0 for False)
"(str)" means that this is a string (but not enclosed in quotes!) 


Basic Program Parameters

#define MAXPOPS   4      // (int) number of populations assumed
#define BURNIN    10000   // (int) length of burnin period
#define NUMREPS   20000   // (int) number of MCMC reps after burnin

Input/Output files

#define INFILE   infile   // (str) name of input data file
#define OUTFILE  outfile  //(str) name of output data file

Data file format

#define NUMINDS    {}    // (int) number of diploid individuals in data file
#define NUMLOCI    {}  // (int) number of loci in data file
#define PLOIDY       1    // (int) ploidy of data
#define MISSING     -9    // (int) value given to missing genotype data
#define ONEROWPERIND 0    // (B) store data for individuals in a single line


#define LABEL     1     // (B) Input file contains individual labels
#define POPDATA   1     // (B) Input file contains a population identifier
#define POPFLAG   0     // (B) Input file contains a flag which says 
                              whether to use popinfo when USEPOPINFO==1
#define LOCDATA   0     // (B) Input file contains a location identifier

#define PHENOTYPE 0     // (B) Input file contains phenotype information
#define EXTRACOLS 0     // (int) Number of additional columns of data 
                             before the genotype data start.

#define MARKERNAMES      0  // (B) data file contains row of marker names
#define RECESSIVEALLELES 0  // (B) data file contains dominant markers (eg AFLPs)
                            // and a row to indicate which alleles are recessive
#define MAPDISTANCES     0  // (B) data file contains row of map distances 
                            // between loci


Advanced data file options

define PHASED           1 // (B) Data are in correct phase (relevant for linkage model only)
#define PHASEINFO        0 // (B) the data for each individual contains a line
                                  indicating phase (linkage model)
#define MARKOVPHASE      0 // (B) the phase info follows a Markov model.
#define NOTAMBIGUOUS  -999 // (int) for use in some analyses of polyploid data

"""


def main():
    """ """
    try:
        infile, outstem = sys.argv[1:]
    except:
        sys.stderr.write("Usage: vcf_to_STRUCTURE.py [in.vcf] [out-file-stem]\n")
        sys.exit(1)
    recs = vcf.Reader(open(infile))
    r = next(recs)
    geno = [ [c.sample, "1", c.gt_nums] for c in r.samples]
    nrecs = 1
    for r in recs:
        # i.e. , is this polymorphic within these samples (not just a fixed
        # non-ref allele
        if r.nucl_diversity > 0:
            for i,gt in enumerate([c.gt_nums for c in r.samples]):
                geno[i].append(gt)
            nrecs += 1
    with open(outstem + ".structure", "w") as out:
        for sample in geno:
            out.write("\t".join(sample) + "\n")

    with open(outstem + ".mparam", "w") as out:
        out.write(params_template.format(len(geno), nrecs))
    print("Write {0} loci for {1} inds to {2}.structure and params file to {2}.mparam".format(nrecs, len(geno), outstem))

if __name__ == "__main__":
    main()
    sys.exit(0)

