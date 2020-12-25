#!/bin/bash -x

for N in 128 256 512
do
  for L in 1 3.14159
  do
    for P in 1 10 20 40
    do
      mpisubmit.pl -p $P  -w 00:15 wave -- $L 0.025 $N 20 POLUS
    done
  done
done
