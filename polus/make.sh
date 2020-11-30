#!/bin/bash -x
#mpixlC -qsmp=omp block.cpp process.cpp function.cpp main.cpp -o main_omp
mpixlC src/*.cpp -o wave
mpixlC src/*.cpp -o wavePi