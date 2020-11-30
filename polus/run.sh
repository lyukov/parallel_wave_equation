#!/bin/bash -x
mpisubmit.pl -p 1  -w 00:15 wave -- 1 0.1 128 20 1 1 1
mpisubmit.pl -p 10 -w 00:15 wave -- 1 0.1 128 20 5 2 1
mpisubmit.pl -p 20 -w 00:15 wave -- 1 0.1 128 20 5 2 2
mpisubmit.pl -p 40 -w 00:15 wave -- 1 0.1 128 20 5 4 2