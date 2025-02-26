#!/bin/bash
#PBS -N map_profile
#PBS -A UCSG0002
#PBS -l walltime=01:00:00
#PBS -q develop
#PBS -j oe
#PBS -l select=1:ncpus=1:mem=4GB

module load linaro-forge
module load mkl

### Run the executable
# cd /glade/work/ctowery/Hackathon-ODE-Miniapp/test
export FORGE_SAMPLER_NUM_SAMPLES=10000
export FORGE_SAMPLER_INTERVAL=2

for i in 1
2
3;
do
cat > user_inputs.nml << EOS
!> Program time integration parameters
&params
    start_time = 0.0,
        !! start time of simulation [s]
    end_time = 2.0,
        !! end time of simulation [s]
    save_name = 'tracers.hst',
        !! file name for full solution outputs
    dt_save = 5e-2,
        !! solution save rate [s]
    nx = 16, 16, 16
        !! 3D domain size
    nflat = 
        !! number of dimensions to flatten into chem ode, choice: {1, 2, 3}
/

!> Initial conditions for carbonate chemistry model
&carbonate_ic
    y_0(1) = 16.0,
        !! Carbon Dioxide, [CO2]... Colin just picked some numbers away from equilibrium
    y_0(2) = 1.67006e03,
        !! Bicarbonate, [HCO3-]
    y_0(3) = 3.14655e02,
        !! Carbonate, [CO32-]
    y_0(4) = 2.96936e02,
        !! Boric Acid, [B(OH)3]
    y_0(5) = 1.18909e02,
        !! Tetrahydroxyborate
    y_0(6) = 20.0
        !! Hydroxide, [OH-]...  Colin just picked some numbers away from equilibrium
    temperature = 25.0,
        !! temperature [degC]
    salinity = 36.9
        !! salinity [PSU]
/
EOS
map --profile --no-mpi --o intel_opt2_flat_1cpu.map ./miniapp_intel_opt2
done
