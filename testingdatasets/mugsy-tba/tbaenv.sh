#!/bin/sh

export MUGSY_INSTALL="/home/user/mugsy"
export MUSCLE_INSTALL="/home/user/muscle"
export TBA_INSTALL="/home/user/multiz-tba"
export PATH=$PATH:$MUGSY_INSTALL:$MUGSY_INSTALL/mapping
export PERL5LIB=$MUGSY_INSTALL/perllibs

TBA_SCRIPT="/home/user/tba"	#directory containing mugsy-tba script
export PATH=$PATH:$TBA_SCRIPT

