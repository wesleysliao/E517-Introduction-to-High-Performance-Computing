#!/bin/bash
#PBS -N mbperf
#PBS -l walltime=00:05:00
#PBS nodes=2:ppn=32
#PBS -q debug_cpu
#PBS -o trace.out

#echo Working directory is $PBS_O_WORKDIR
#cd $PBS_O_WORKDIR

RUN_COMMAND=aprun
RUN_COMMAND=/usr/lib64/openmpi/bin/mpirun
MAX_PS_POWER2=6
MAX_PS_POWER2=2

SS_PIXELS=32768
SS_PIXELS=512
WS_PIXELS=8192
WS_PIXELS=128
PERF_PIXELS=$SS_PIXELS

# 
# Strong Scaling
#
echo STONG SCALING TEST

for i in $(seq 0 $MAX_PS_POWER2); do
    processors=$((2 ** i))
    echo computing mandelbrot of $SS_PIXELS by $SS_PIXELS \($((SS_PIXELS * SS_PIXELS))\) on $processors processors

    tout=$(time $RUN_COMMAND -n $processors ./mb $SS_PIXELS $SS_PIXELS)
    echo $tout
done

#
# Weak Scaling
#
echo WEAK SCALING TEST

for i in $(seq 0 $MAX_PS_POWER2); do
    processors=$((2 ** i))
    echo computing mandelbrot of $WS_PIXELS by $((processors * WS_PIXELS)) \($((WS_PIXELS * processors * WS_PIXELS))\) on $processors processors

    $RUN_COMMAND -n $processors ./mb $WS_PIXELS $((processors * WS_PIXELS))
done

#
# Perf results
#
echo PERF TEST
echo $((2 ** $MAX_PS_POWER2 - 1))
$RUN_COMMAND -n $((2 ** $MAX_PS_POWER2 - 1)) ./mb $PERF_PIXELS $PERF_PIXELS : -n 1 perf record ./mb $PERF_PIXELS $PERF_PIXELS
