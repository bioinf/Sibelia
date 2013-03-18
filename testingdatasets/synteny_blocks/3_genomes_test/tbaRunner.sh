#!/bin/bash

TREE="(Xuzhou21,(UTI89,W))"

if (( $# < 2 )); then
	echo "args, please!"
	exit
fi

while read LINE; do
	ARGS=( $LINE )
	DIRNAME="tba_${ARGS[0]}_${ARGS[1]}"
	echo $DIRNAME:
	if [[ ! -d $DIRNAME ]]; then
		mkdir $DIRNAME
	fi
	./mugsy-tba --tba --tree $TREE --directory `pwd`/$DIRNAME --distance ${ARGS[0]} --minlength ${ARGS[1]} `cat $2`
	echo ${ARGS[1]} > "$DIRNAME/blocksize.txt"
done < $1
