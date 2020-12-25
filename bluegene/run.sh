#!/bin/bash -x


for N in 128 256 512
do
  for L in 1 3.14159
  do
    mpisubmit.bg -n 1   -w 00:15:00 -m smp wave -- $L 0.025 $N 20 BLUEGENE
    mpisubmit.bg -n 64  -w 00:15:00 -m smp wave -- $L 0.025 $N 20 BLUEGENE
    mpisubmit.bg -n 128 -w 00:15:00 -m smp wave -- $L 0.025 $N 20 BLUEGENE
    mpisubmit.bg -n 256 -w 00:10:00 -m smp wave -- $L 0.025 $N 20 BLUEGENE
    mpisubmit.bg -n 512 -w 00:05:00 -m smp wave -- $L 0.025 $N 20 BLUEGENE
  done
done

for N in 128 256 512
do
  for L in 1 3.14159
  do
    mpisubmit.bg -n 1   -w 00:15:00 -m smp -e "OMP_NUM_THREADS=2" wave_omp -- $L 0.025 $N 20 BLUEGENE_OMP
    mpisubmit.bg -n 64  -w 00:15:00 -m smp -e "OMP_NUM_THREADS=2" wave_omp -- $L 0.025 $N 20 BLUEGENE_OMP
    mpisubmit.bg -n 128 -w 00:15:00 -m smp -e "OMP_NUM_THREADS=2" wave_omp -- $L 0.025 $N 20 BLUEGENE_OMP
    mpisubmit.bg -n 256 -w 00:10:00 -m smp -e "OMP_NUM_THREADS=2" wave_omp -- $L 0.025 $N 20 BLUEGENE_OMP
    mpisubmit.bg -n 512 -w 00:05:00 -m smp -e "OMP_NUM_THREADS=2" wave_omp -- $L 0.025 $N 20 BLUEGENE_OMP
  done
done
