module load SpectrumMPI
module load OpenMPI
nvcc -ccbin mpicxx -std=c++11 -O3 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_60,code=sm_60 \
     src/*.cu src/*.cpp -o cuda_wave
mpixlC -std=c++11 -qsmp=omp src/*.cpp -o omp_wave
mpixlC -std=c++11 src/*.cpp -o simple_wave