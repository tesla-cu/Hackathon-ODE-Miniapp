#!/bin/bash
#PBS -N acc_test
#PBS -A UCSG0002
#PBS -l walltime=01:00:00
#PBS -q main
#PBS -j oe
#PBS -l select=1:ncpus=16:ompthreads=16

module load nvhpc

### Run the executable
export ACC_NUM_CORES=16
# cd /glade/work/ctowery/Hackathon-ODE-Miniapp/test
./miniapp.exe
