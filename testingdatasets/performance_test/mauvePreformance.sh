#!/bin/bash 

#Command line used to run mauve on big set of genomes
#Change paths if needed

MAUVE_EXEC="./progressiveMauve"
INPUT_PATH="testgenomes"

$MAUVE_EXEC --output=mauve.xmfa `find $INPUT_PATH -type f`
