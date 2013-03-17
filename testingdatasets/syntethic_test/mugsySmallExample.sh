#!/bin/bash 

#Command line used to run mugsy on syntenic exmple (see results section)
#This file should be placed in mugsy`s installation directory.
#Don`t forget to prepare mugsy (source mugsyenv.sh etc.)

MUGSY_EXEC="./mugsy"
INPUT_PATH="testgenomes"

$MUGSY_EXEC --directory=. $INPUT_PATH/smallTest1.fasta $INPUT_PATH/smallTest2.fasta
