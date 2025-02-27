#!/bin/bash
#PBS -N acc_test
#PBS -A UCSG0002
#PBS -l walltime=01:00:00
#PBS -q tutorial
#PBS -j oe
#PBS -l select=1:ncpus=64:mpiprocs=1:ngpus=1

module load nvhpc

### Run the executable
# cd /glade/work/ctowery/Hackathon-ODE-Miniapp/test
./miniapp.exe
