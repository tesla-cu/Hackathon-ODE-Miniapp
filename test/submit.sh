#!/bin/bash
#PBS -N map_profile
#PBS -A UCSG0002
#PBS -l walltime=01:00:00
#PBS -q develop
#PBS -j oe
#PBS -l select=1:ncpus=1:mem=4GB

module load mkl
# module load cce/17
# module load linaro-forge

export TMPDIR=$SCRATCH/temp
mkdir -p $TMPDIR

### Run the executable
# cd /glade/work/ctowery/Hackathon-ODE-Miniapp/test
export FORGE_SAMPLER_NUM_SAMPLES=10000
export FORGE_SAMPLER_INTERVAL=2
map --connect --no-mpi ./miniapp.exe
