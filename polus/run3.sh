#!/bin/bash -x
mpisubmit.pl -p 1  -w 00:15 wavePi -- 3.14 0.1 129 20 1 1 1
mpisubmit.pl -p 10 -w 00:15 wavePi -- 3.14 0.1 129 20 5 2 1
mpisubmit.pl -p 20 -w 00:15 wavePi -- 3.14 0.1 129 20 5 2 2
mpisubmit.pl -p 40 -w 00:15 wavePi -- 3.14 0.1 129 20 5 4 2

mpisubmit.pl -p 1  -w 00:15 wavePi -- 3.14 0.05 257 20 1 1 1
mpisubmit.pl -p 10 -w 00:15 wavePi -- 3.14 0.05 257 20 5 2 1
mpisubmit.pl -p 20 -w 00:15 wavePi -- 3.14 0.05 257 20 5 2 2
mpisubmit.pl -p 40 -w 00:15 wavePi -- 3.14 0.05 257 20 5 4 2

mpisubmit.pl -p 1 -w 00:15  wavePi -- 3.14 0.025 513 20 1 1 1
mpisubmit.pl -p 10 -w 00:15 wavePi -- 3.14 0.025 513 20 5 2 1
mpisubmit.pl -p 20 -w 00:15 wavePi -- 3.14 0.025 513 20 5 2 2
mpisubmit.pl -p 40 -w 00:15 wavePi -- 3.14 0.025 513 20 5 4 2