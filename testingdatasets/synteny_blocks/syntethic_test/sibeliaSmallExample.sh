#!/bin/bash

#Command line used to run Sibelia on syntenic exmple (see results section)

SIB_EXEC="./Sibelia"
INPUT_PATH="./testgenomes"

$SIB_EXEC -s fine $INPUT_PATH/smallTest1.fasta $INPUT_PATH/smallTest2.fasta
