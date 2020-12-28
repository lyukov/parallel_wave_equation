module load SpectrumMPI
module load OpenMPI

for L in 1.0 3.14
do
    for N in 128 256 512
    do
        for P in 1 2 3 4
        do
        mpisubmit.pl --stdout log/simple/$L"_"$N"_"$P.out --stderr log/simple/$L"_"$N"_"$P.err -p $P -w 00:15 ./simple_wave -- $L 0.025 $N 20 simple
        done
    done
done