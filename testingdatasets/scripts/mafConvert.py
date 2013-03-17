#!/usr/bin/env python

import sys
import os

sys.path.insert(1, os.path.join(os.path.dirname(__file__),"alignio-maf"))
from Bio import AlignIO

if len(sys.argv) <= 3:
	print "Usage: mafConvert maf_file min_block_size strains_file"
	sys.exit(1)

seqNames = {}
seqCount = 1
output = ""
MIN_BLOCK_SIZE = int(sys.argv[2])
FILTER = True

#read strains file 
for line in open(sys.argv[3], "r"):
	#print line
	line = line.strip()
	vals = line.split(" ")
	seqNames[vals[1]] = vals[0]
#print seqNames

counter = 1
for block in AlignIO.parse(sys.argv[1], "maf"):
	if block[0].annotations["size"] < MIN_BLOCK_SIZE:
		continue

	if FILTER and len(block) == 1:	#multiplicity 1
		continue

	output += "#{0}\n".format(str(counter))
	for seq in block:
		name = seq.name.split('.')[0]
		#print name, seqNames
		assert(name in seqNames)
		seqId = seqNames[name]
				
		output += str(seqId) + "\t" + seq.annotations["strand"] + "\t"
		output += str(seq.annotations["start"]) + "\t" + str(seq.annotations["size"]) + "\n"
	output += "\n"
	counter += 1

for key in sorted(seqNames, key=seqNames.get):
	print str(seqNames[key]) + "\t" + key
print ""

print output
