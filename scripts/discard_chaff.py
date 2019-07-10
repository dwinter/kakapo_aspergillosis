#!/usr/bin/python

#
# Take a genome possible including many small contigs and filter out 
# those smaller than a given size.
#
#  Usage: 
#    discard_chaff.py [scaffolds.fa] [min-len]
# 
# contigs with length greater than or equal to the min-len are printed to stdout

import sys
from Bio import SeqIO


def skip_short_contigs(recs, n=500):
    """GENERATOR FUNCTION
    
    Iterates over a set of SeqRecord objects and yields only those
    with length > n.
    """
    for r in recs:
        if len(r) >= n:
            yield(r)

if __name__ == "__main__":
    try:
        contigs, n = sys.argv[1:]
    except ValueError:
        print("Usage: discard_chaff.py [scaffolds.fa] [min-len]")
        sys.exit(1)
    wheat = skip_short_contigs(SeqIO.parse(contigs, "fasta"), int(n))
    SeqIO.write(wheat, sys.stdout, "fasta")
    sys.exit(0)
