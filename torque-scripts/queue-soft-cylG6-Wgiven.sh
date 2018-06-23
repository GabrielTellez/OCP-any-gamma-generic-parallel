#!/bin/sh
W=$1
moreargs=$2
for n in 2 3 4 5 6 7 8 9 10 11 12 13 14 ; do
./ocp -b 100000 -t 48 --model 4 -i ../binaryfiles/gamma3_n${n}.bin -p $W --compute-coefs $moreargs >> res/soft_cyl_resultsG6_W${W}.txt ;
done 


