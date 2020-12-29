module load SpectrumMPI
module load OpenMPI

for L in 1.0 3.14
do
    for N in 128 256 512
    do
        for P in 1 2 3 4
        do
        echo "source /polusfs/setenv/setup.SMPI" > $L"_"$N"_"$P".lsf"
        echo "#BSUB -n "$P >> $L"_"$N"_"$P".lsf"
        echo "#BSUB -W 00:15" >> $L"_"$N"_"$P".lsf"
        echo "#BSUB -o log/omp/"$L"_"$N"_"$P".out" >> $L"_"$N"_"$P".lsf"
        echo "#BSUB -e log/omp/"$L"_"$N"_"$P".err" >> $L"_"$N"_"$P".lsf"
        echo "#BSUB -R \"span[ptile=1]\"" >> $L"_"$N"_"$P".lsf"
        echo "OMP_NUM_THREADS=8 mpiexec ./omp_wave "$L" 0.025 "$N" 40 omp" >> $L"_"$N"_"$P".lsf"
        bsub < $L"_"$N"_"$P".lsf"
        rm $L"_"$N"_"$P".lsf"
        done
    done
done