#!/bin/bash

GENES=`ls $1`
ARGS=( $GENES )

STRING="${ARGS[0]%.fasta}"

for I in `seq 1 $(( ${#ARGS[@]} - 1 ))`; do
	FILE=${ARGS[I]}
	#FILE=${FILE%.fasta}
	STRING="(${FILE%.fasta},${STRING})"
done

echo $STRING
