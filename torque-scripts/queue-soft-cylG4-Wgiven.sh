#!/bin/sh
W=$1
moreargs=$2
script=run-soft-cylG4-W$W
cat > $script <<EOF
#!/bin/sh
#
#PBS -l nodes=1:ppn=8
#PBS -q fisica
#
WORKDIR=/fisica/gtellez/gamma4-6/any-gamma/generic-parallel
cd \$WORKDIR
for n in 2 3 4 5 6 7 8 9 10 11 12 13 14 ; do
./ocp -b 100000 -t 48 --model 4 -i ../binaryfiles/bosons-n\${n}.bin -p $W --compute-coefs $moreargs >> res/soft_cyl_resultsG4_W${W}.txt ;
done 
EOF
cat $script
#qsub $script


