#!/usr/bin/env python

import sys
from bitarray import bitarray
from itertools import repeat
import os

if len(sys.argv) <= 2:
	print "Usage: boundaryDiscord.py first_file second_file"
	sys.exit(1)

def retriveBlocks(filename):
	inFile = open(filename)
	numSeq = 0
	while True:
		line = inFile.readline()
		line = line.strip()
		if len(line) == 0:
			break
		numSeq += 1

	blocks = []
	for i in xrange(numSeq):
		blocks.append([])

	while True:
		line = inFile.readline()
		line = line.strip()
		if len(line) == 0:
			break
		assert(line[0] == '#') #block id

		while True:
			seq = inFile.readline()
			seq = seq.strip()
			if (len(seq) == 0):
				break
				
			vals = seq.split('\t')
			seqid = int(vals[0])
			start = int(vals[2])
			length = int(vals[3])
			
			blocks[seqid - 1].append((start, length))
	for block in blocks:
		block.sort(key = lambda b : b[0])
	return blocks

trueBlocks = retriveBlocks(sys.argv[1])
checkBlocks = retriveBlocks(sys.argv[2])


for i in xrange(len(trueBlocks)):
	sumScore = 0.0
	sumNormScore = 0.0
	counter = 0

	for block1 in trueBlocks[i]:
		bestScore = float("-inf")
		num = 0
		for block2 in checkBlocks[i]:
			if not block2:
				continue
			
			if block2[0] > block1[0] + block1[1]:
				break

			left = max(block1[0], block2[0])
			right = min(block1[0] + block1[1], block2[0] + block2[1])
			overlap = right - left if left <= right else 0
			if overlap > block1[1] / 2 and overlap > block2[1] / 2:
				#print str(block1), "overlaps", str(block2), overlap
				#score = abs(block1[0] - block2[0]) + abs(block1[0] + block1[1] - block2[0] - block2[1])
				score = overlap
				num += 1
				if score > bestScore:
					#bestOverlap = overlap
					bestScore = score
					bestBlock = block2
		print num
		#assert(num <=1 )

		if bestScore != float("-inf"):
			left = min(block1[0], bestBlock[0])
			right = max(block1[0] + block1[1], bestBlock[0] + bestBlock[1])
			sys.stderr.write(str(block1) + " overlaps " + str(bestBlock) + " with score " + str(bestScore))
			sys.stderr.write(" " + str( float(bestScore) / (right - left) ) + "\n")
			sumScore += bestScore ** 2
			sumNormScore += (float(bestScore) / (right - left)) ** 2
			counter += 1
			#checkBlocks[checkBlocks.index(block2)] = None

	if counter > 0:
		sys.stdout.write(str(i + 1) + ": ")
		sys.stdout.write("Overlapped: " + str(counter) + " Mean score: " + str((sumScore / counter) ** 0.5))
		sys.stdout.write(" Normalized: " + str((sumNormScore / counter) ** 0.5) + "\n")
