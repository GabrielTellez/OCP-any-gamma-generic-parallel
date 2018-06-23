#!/bin/bash
Gamma=$1
beg=$2
en=$3
moreargs=$4
for (( c=$beg; c<=$en; c++ ))
do
    for (( d=0; d<=9; d++ ))
    do
	./queue-soft-cylG${Gamma}-Wgiven.sh $c.$d $moreargs
    done
done
