module load SpectrumMPI
module load OpenMPI

for L in 1.0 3.14
do
    for N in 128 256 512
    do
        for P in 1 2 3 4
        do
        export OMP_NUM_THREADS=8
        bsub -q normal -W 00:15 -x -R "span[ptile=1]" -n $P \
             -oo log/omp/${L}_${N}_${P}.out \
             -eo log/omp/${L}_${N}_${P}.err \
             mpirun -n ${P} ./omp_wave $L 0.025 $N 40 omp
        done
    done
done