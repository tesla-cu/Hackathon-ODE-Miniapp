#!/bin/bash
#PBS -N nvtest
#PBS -A UCSG0002
#PBS -l walltime=01:00:00
#PBS -q develop
#PBS -j oe
#PBS -l select=1:ncpus=1:mem=8GB

module load nvhpc

### Run the executable
# cd /glade/work/ctowery/Hackathon-ODE-Miniapp/test
./miniapp.exe
