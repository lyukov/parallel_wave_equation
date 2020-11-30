#!/bin/bash -x
#mpixlC -qsmp=omp block.cpp process.cpp function.cpp main.cpp -o main_omp
mpixlC *.cpp -o wave
mpixlC *.cpp -o wavePi