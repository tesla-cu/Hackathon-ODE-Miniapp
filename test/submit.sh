#!/bin/bash
#PBS -N miniapp_test
#PBS -A UCSG0002
#PBS -l walltime=01:00:00
#PBS -q develop
#PBS -j oe
#PBS -l select=1:ncpus=1:mem=4GB

module load linaro-forge

export TMPDIR=$SCRATCH/temp
mkdir -p $TMPDIR

### Run the executable
cd /glade/work/ctowery/Hackathon-ODE-Miniapp/test
map --connect ./miniapp.exe
