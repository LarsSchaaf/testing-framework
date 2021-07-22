#!/bin/bash

#PBS -l select=1:ncpus=40:mpiprocs=40:ompthreads=1:mem=170GB
#PBS -l walltime=12:10:00
#PBS -q computeq
cd ${PBS_O_WORKDIR}/

module purge
module load intel/compiler/2017.5
module load rt_intel_impi/2017.5
module load intel/python3.6/2019.1
module load quantum_espresso/6.5/2020-02-13_09-05-31
mkdir scratch
export SCRATCH=${PBS_O_WORKDIR}/scratch/
export PYTHONPATH=$PYTHONPATH:/gpfs/nobackup/projects/quantumchemistry/rom_cq/SOAPML/venv/SOAPML/lib/python3.6/site-packages/



python run.py ./phonopy_R3c_traj.xyz >>out
#python run.py >>out
