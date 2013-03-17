#!/bin/bash

#Command line used to run preformance test for Sibelia

SIB_EXEC="./Sibelia"
INPUT_PATH="./testgenomes"

$SIB_EXEC -s fine `find $INPUT_PATH -type f`
