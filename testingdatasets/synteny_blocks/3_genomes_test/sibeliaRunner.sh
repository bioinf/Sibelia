#!/bin/bash

if (( $# < 2 )); then
	echo "args, please!"
	exit
fi

while read LINE; do
	ARGS=( $LINE )
	if (( ${#ARGS[@]} < 3 )); then
		DIRNAME="sibelia_${ARGS[0]}_${ARGS[1]}"
		echo $DIRNAME:
		if [[ ! -d $DIRNAME ]]; then
			mkdir $DIRNAME
		fi

		./Sibelia -s ${ARGS[0]} -m ${ARGS[1]} -o $DIRNAME `cat $2`
	else
		DIRNAME="sibelia_${ARGS[0]}_${ARGS[1]}_${ARGS[2]}"
		echo $DIRNAME:
		if [[ ! -d $DIRNAME ]]; then
			mkdir $DIRNAME
		fi
		echo -e "2\n30 150" > stagefile.txt
		echo ${ARGS[0]} ${ARGS[1]} >> stagefile.txt
		./Sibelia -k stagefile.txt -m ${ARGS[2]} -o $DIRNAME `cat $2`  
	fi
	done < $1
