#!/bin/bash -x

for N in 128 256 512
do
  mpisubmit.bg -n 1   -w 00:15:00 -m smp wavePi -- 3.14159 0.05 $N 20 1 1 1
  mpisubmit.bg -n 64  -w 00:15:00 -m smp wavePi -- 3.14159 0.05 $N 20 4 4 4
  mpisubmit.bg -n 128 -w 00:15:00 -m smp wavePi -- 3.14159 0.05 $N 20 8 4 4
  mpisubmit.bg -n 256 -w 00:10:00 -m smp wavePi -- 3.14159 0.05 $N 20 8 8 4
done
