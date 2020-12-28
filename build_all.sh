module load SpectrumMPI
module load OpenMPI
nvcc -ccbin mpixlC -O3 -std=c++11 src/*.cu src/*.cpp -o cuda_wave
mpixlC -std=c++11 -qsmp=omp src/*.cpp -o omp_wave
mpixlC -std=c++11 src/*.cpp -o simple_wave