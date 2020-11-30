#!/bin/bash -x

for N in 128 256 512
do
  mpisubmit.pl -p 1  -w 00:15 wave -- 1 0.05 $N 20 1 1 1
  mpisubmit.pl -p 10 -w 00:15 wave -- 1 0.05 $N 20 5 2 1
  mpisubmit.pl -p 20 -w 00:15 wave -- 1 0.05 $N 20 5 2 2
  mpisubmit.pl -p 40 -w 00:15 wave -- 1 0.05 $N 20 5 4 2
done

for N in 128 256 512
do
  mpisubmit.pl -p 1  -w 00:15 wavePi -- 3.14159 0.1 $N 20 1 1 1
  mpisubmit.pl -p 10 -w 00:15 wavePi -- 3.14159 0.1 $N 20 5 2 1
  mpisubmit.pl -p 20 -w 00:15 wavePi -- 3.14159 0.1 $N 20 5 2 2
  mpisubmit.pl -p 40 -w 00:15 wavePi -- 3.14159 0.1 $N 20 5 4 2
done