#!/bin/bash
mpic++ -std=c++98 -Xpreprocessor -fopenmp -O3 src/*.cpp -o bin/wave -lomp