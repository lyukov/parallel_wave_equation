module load SpectrumMPI
module load OpenMPI
nvcc -ccbin mpixlC -O3 -std=c++11 src/*.cu src/*.cpp -o cuda_wave