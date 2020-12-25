nvcc -ccbin mpicxx -std=c++11 -O3 -Xcompiler "-g -pg" \
     -gencode arch=compute_52,code=sm_52 -gencode arch=compute_60,code=sm_60 src/*.cu src/*.cpp -o cuda_wave && \
mpisubmit.pl -p 1 -g -- ./cuda_wave 1 0.025 512 20 1 1 1 cuda && \
gprof cuda_wave gmon.out > gprof.stats