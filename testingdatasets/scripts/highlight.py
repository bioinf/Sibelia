#!/usr/bin/env python

import sys
import os
import re

if len(sys.argv) <= 1:
	print "Usage: highlight.py file"
	sys.exit(1)

for line in open(sys.argv[1], "r"):
	m = re.match("([0-9]+)\t[+-][0-9]+\t([0-9]+)\t([0-9]+)", line)
	if m:
		print "seq" + m.group(1), m.group(2), str(int(m.group(3)) + int(m.group(2)))
