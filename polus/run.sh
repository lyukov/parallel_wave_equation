#!/bin/bash -x

for N in 128 256 512
do
  for L in 1 3.14159
  do
    mpisubmit.pl -p 1  -w 00:15 wave -- $L 0.1 $N 20 1 1 1 POLUS
    mpisubmit.pl -p 10 -w 00:15 wave -- $L 0.1 $N 20 5 2 1 POLUS
    mpisubmit.pl -p 20 -w 00:15 wave -- $L 0.1 $N 20 5 2 2 POLUS
    mpisubmit.pl -p 40 -w 00:15 wave -- $L 0.1 $N 20 5 4 2 POLUS
  done
done
