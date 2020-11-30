#!/bin/bash -x
mpixlcxx_r -qsmp=omp src/*.cpp -o wave_omp
mpixlcxx_r           src/*.cpp -o wave
mpixlcxx_r -qsmp=omp src/*.cpp -o wavePi_omp
mpixlcxx_r           src/*.cpp -o wavePi
