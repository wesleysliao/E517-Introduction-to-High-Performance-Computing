#!/bin/bash
#PBS -N hd 
#PBS -l walltime=00:05:00
#PBS nodes=2:ppn=32
#PBS -q debug_cpu
#PBS -o trace.out

#RUN CONDITIONS
#echo Working directory is $PBS_O_WORKDIR
#cd $PBS_O_WORKDIR
RUN_COMMAND=aprun
MAX_PS_POWER2=6

#TEST CONDITIONS
RUN_COMMAND=/usr/lib64/openmpi/bin/mpirun
MAX_PS_POWER2=2

for i in $(seq 0 $MAX_PS_POWER2); do
    processors=$((2 ** i))
    echo heat distribution with $processors processors

    $RUN_COMMAND -n $processors ./hd
done
