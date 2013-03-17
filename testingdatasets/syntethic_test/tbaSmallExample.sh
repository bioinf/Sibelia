#!/bin/bash

#Command line used to run tba on syntenic exmple (see results section)
#This file should be placed in tba`s installation directory.
#Make sure, that mugsy-tba is installed properly (see README for instructions)

TBA_EXEC="./mugsy-tba"
INPUT_PATH="testgenomes"

$TBA_EXEC --tba --tree "(smallTest1,smallTest2)" --directory . $INPUT_PATH/smallTest1.fasta \
	$INPUT_PATH/smallTest2.fasta
