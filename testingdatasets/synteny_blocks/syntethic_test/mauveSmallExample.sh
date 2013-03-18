#!/bin/bash 

#Command line used to run mauve on syntenic exmple (see results section)
#Change paths if needed

MAUVE_EXEC="./progressiveMauve"
INPUT_PATH="testgenomes"

$MAUVE_EXEC --output=mauve.xmfa $INPUT_PATH/smallTest1.fasta $INPUT_PATH/smallTest2.fasta
