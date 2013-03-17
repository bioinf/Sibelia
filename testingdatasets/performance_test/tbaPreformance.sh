#!/bin/bash

#Command line used to run tba on big set of genomes
#This file should be placed in tba`s installation directory.
#Make sure, that mugsy-tba is installed properly (see README for instructions)

TBA_EXEC="./mugsy-tba"
INPUT_PATH="testgenomes"

$TBA_EXEC --tba --tree `bash guidetree.sh $INPUT_PATH` --directory . `find $INPUT_PATH -type f`