#!/bin/bash

if (( $# < 2 )); then
	echo "args, please!"
	exit
fi

while read LINE; do
	ARGS=( $LINE )
	DIRNAME="mauve_${ARGS[0]}"
	echo $DIRNAME:
	if [[ ! -d $DIRNAME ]]; then
		mkdir $DIRNAME
	fi
	./progressiveMauve --output=`pwd`/$DIRNAME --weight=${ARGS[0]} `cat $2`
	echo ${ARGS[0]} > "$DIRNAME/blocksize.txt"
done < $1
