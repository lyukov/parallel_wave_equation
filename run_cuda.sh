module load SpectrumMPI
#module load OpenMPI

for L in 1.0 3.14
do
    for N in 128 256 512
    do
        for P in 1 2 3 4
        do
        bsub -q normal -W 00:10 -x -R "span[ptile=1]" -gpu "num=1" -n $P \
             -oo log/cuda/$L"_"$N"_"$P.out -eo log/cuda/$L"_"$N"_"$P.err \
             mpirun -n $P ./cuda_wave $L 0.025 $N 40 cuda
        done
    done
done