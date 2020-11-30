#!/bin/bash -x
#mpixlC -qsmp=omp block.cpp process.cpp function.cpp main.cpp -o main_omp
module load SpectrumMPI
mpixlC src/*.cpp -o wave