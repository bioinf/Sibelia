#!/bin/bash 

#Command line used to run mugsy on big set of genomes
#This file should be placed in mugsy`s installation directory.
#Don`t forget to prepare mugsy (source mugsyenv.sh etc.)

MUGSY_EXEC="./mugsy"
INPUT_PATH="testgenomes"

$MUGSY_EXEC --directory=. `find $INPUT_PATH -type f`
