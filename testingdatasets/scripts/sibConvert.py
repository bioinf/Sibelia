#!/usr/bin/env python

import sys

if len(sys.argv) <= 1:
	print "Usage: subConvert.py sib_file"
	sys.exit(1)

fileIn = open(sys.argv[1])

fileIn.readline() #header
while True:
	line = fileIn.readline().strip()
	if line[0] == '-':
		break
	vals = line.split('\t')
	print vals[0] + '\t' + vals[2]

print ""

while True:
	line = fileIn.readline().strip()
	if not line:
		break
	print line.split(' ')[1]
	fileIn.readline()
	while True:
		line = fileIn.readline().strip()
		if line[0] == '-':
			print ""
			break
		vals = line.split('\t')
		print vals[0] + "\t" + vals[1] + "1\t" + vals[2] + "\t" + vals[4]
