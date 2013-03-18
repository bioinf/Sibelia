#!/usr/bin/env python

import sys
from bitarray import bitarray
from itertools import repeat
import os

if len(sys.argv) <= 2:
	print "Usage: f1stat.py first_file second_file"
	sys.exit(1)
	
def retriveCoverage(filename):
	inFile = open(filename)
	numSeq = 0
	while True:
		line = inFile.readline()
		line = line.strip()
		if len(line) == 0:
			break
		numSeq += 1

	totalCover = []
	for i in xrange(numSeq):
		totalCover.append(bitarray())

	bcounter = 0

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
			blockid = int(vals[0]) - 1
			start = int(vals[2])
			length = int(vals[3])

			for j in xrange(start, start + length + 1):
				l = totalCover[blockid].length()
				if j >= l:
					totalCover[blockid].extend(repeat(False, j - l + 1))
				totalCover[blockid][j] = True
		sys.stderr.write(filename + ": block " + str(bcounter) + "\n")
		bcounter += 1
	
	return totalCover

trueCover = retriveCoverage(sys.argv[1])
checkCover = retriveCoverage(sys.argv[2])

for i in xrange(len(trueCover)):
	falsePositive = 0
	falseNegative = 0
	truePositive = 0
	trueNegative = 0

	for j in xrange(len(trueCover[i])):
		if len(trueCover[i]) > len(checkCover[i]):
			checkCover[j].extend(repeat(False, len(trueCover[i]) - len(checkCover[i]) + 1))

		if trueCover[i][j]:
			if checkCover[i][j]:
				truePositive += 1
			else:
				falsePositive += 1
		else:
			if checkCover[i][j]:
				falseNegative += 1
			else:
				trueNegative += 1
	sys.stdout.write(str(i + 1) + ":\n")
	sys.stdout.write("tp: {0}, tn: {1}, fp: {2}, fn: {3}\n".format(truePositive, trueNegative, 
																	falsePositive, falseNegative))
	precision = float(truePositive) / (truePositive + falsePositive)
	recall = float(truePositive) / (truePositive + falseNegative)
	f = 2 * precision * recall / (precision + recall)
	sys.stdout.write("P = {0}, R = {1}, F = {2}\n".format(precision, recall, f ))
