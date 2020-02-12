#!/bin/bash
#PBS -N ass41 
#PBS -l walltime=00:05:00
#PBS nodes=2:ppn=32
#PBS -q debug_cpu
#PBS -o trace.out

#RUN CONDITIONS
#echo Working directory is $PBS_O_WORKDIR
#cd $PBS_O_WORKDIR
RUN_COMMAND=aprun
MAX_PS_POWER2=6
PI_DIVISIONS=1000000000
INT_DIVISIONS=1000000000

#TEST CONDITIONS
RUN_COMMAND=/usr/lib64/openmpi/bin/mpirun
MAX_PS_POWER2=2
PI_DIVISIONS=10000000
INT_DIVISIONS=10000000

echo CALCULATION OF PI STRONG SCALING TEST 

for i in $(seq 0 $MAX_PS_POWER2); do
    processors=$((2 ** i))
    echo computing pi with $processors processors

    tout=$(time $RUN_COMMAND -n $processors ./p1 $PI_DIVISIONS)
    echo $tout
done

echo NUMERICAL INTEGRATION STRONG SCALING TEST 

for i in $(seq 0 $MAX_PS_POWER2); do
    processors=$((2 ** i))
    echo computing integral with $processors processors

    tout=$(time $RUN_COMMAND -n $processors ./p2 $INT_DIVISIONS)
    echo $tout
done
