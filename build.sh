module load SpectrumMPI
module load OpenMPI
nvcc -ccbin mpixlC -O3 \
     src/*.cu src/*.cpp -o cuda_wave