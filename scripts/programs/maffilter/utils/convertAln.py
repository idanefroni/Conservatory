#! /usr/bin/python

# Convert an alignment into a single-block MAF file.
import sys
from Bio import AlignIO
alignment = AlignIO.read(open(sys.argv[1] + ".phy"), "phylip")
print "Alignment length %i" % alignment.get_alignment_length()
# Convert to maf:
maf = open(sys.argv[1] + ".maf", 'w')
maf.write("##maf version=1\n")
maf.write("#Generated from a single alignment\n\n")
maf.write("a\n")
for record in alignment :
    maf.write("s seq%s.1\t0\t%d\t+\t%d\t%s\n" % (record.id, len(record.seq), len(record.seq), record.seq))
maf.close()
