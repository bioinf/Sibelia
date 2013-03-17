#!/usr/bin/env python

import sys
from bitarray import bitarray
import os

if len(sys.argv) <= 1:
	print "Usage: coverage.py file"
	sys.exit(1)

inFile = open(sys.argv[1], "r")
GENES_PATH = "../comparison/genomes/"

#get approx sequences lengths
seqSize = []
while True:
	line = inFile.readline()
	line = line.strip()
	if len(line) == 0:
		break
	vals = line.split('\t')
	filename = vals[1].split('.')[0] + ".fasta"
	seqSize.append(os.path.getsize(GENES_PATH + filename))
numSeq = len(seqSize)
coverTable = [numSeq * [0] for _ in xrange(numSeq)]
seqCount = [numSeq * [0]]

totalCover = []
for seqlen in seqSize:
	#totalCover.append([False] * seqlen)
	totalCover.append(bitarray(seqlen))
	totalCover[-1].setall(False)

#covered = [numSeq * False]
bcounter = 0

while True:
	line = inFile.readline()
	line = line.strip()
	if len(line) == 0:
		break
	assert(line[0] == '#') #block id
	ids = []
	start = []
	length = []
	count = 0
	while True:
		seq = inFile.readline()
		seq = seq.strip()
		#end of block
		if (len(seq) == 0):
			break
		count += 1
		vals = seq.split('\t')
		ids.append(int(vals[0]) - 1)
		start.append(int(vals[2]))
		length.append(int(vals[3]))
	
	for i in xrange(count):
		for j in xrange(start[i], start[i] + length[i] + 1):
			if j >= seqSize[ids[i]]:
				print "alarm!", j
				continue
			totalCover[ids[i]][j] = True
	print "block", bcounter
	bcounter += 1

for i in xrange(numSeq):
	print i + 1, ":",
	covered = 0
	for b in totalCover[i]:
		if b:
			covered += 1
	print float(covered) / seqSize[i]
