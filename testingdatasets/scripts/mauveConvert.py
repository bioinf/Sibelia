#!/usr/bin/env python

import sys
import os
import re

sys.path.insert(1, os.path.join(os.path.dirname(__file__),"alignio-maf"))
from Bio import AlignIO

if len(sys.argv) <= 3:
	print "Usage: mauveConvert.py xmfa_file min_block_size strains_file"
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
output += "\n#1\n"
mult = 0
blocks = ""
for line in open(sys.argv[1], "r"):
	if line[0] == '=':
		if FILTER and mult > 1:
			output += blocks
			counter += 1
		blocks = ""
		mult = 0
		if counter != 1:
			blocks += "\n#%d\n" % counter
		mult = 0
	if line[0] != '>':
		continue
	m = re.match("> ([0-9]*):([0-9]*)-([0-9]*) ([+-]) (.*)", line)
	seqid = m.group(1)
	start = m.group(2)
	end = m.group(3)
	strand = m.group(4)
	seqname = m.group(5)
	#strain = re.match(".*/(.*).fasta", seqname)
	#print strain.group(1)
	#seqid = seqNames[strain.group(1)]
	length = int(end) - int(start)

	if length < MIN_BLOCK_SIZE:
		continue

	mult += 1
	blocks += seqid + "\t" + strand + "1\t" + start + "\t" + str(length) + "\n"
	#if not seqname in seqNames:
	#	seqNames[seqname] = seqid
		


for key in sorted(seqNames, key = seqNames.get):
	print str(seqNames[key]) + "\t" + key

print output
