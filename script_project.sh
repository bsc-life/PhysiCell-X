#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH -t 08:00:00
#SBATCH -o output-%j
#SBATCH -e error-%j
#SBATCH --exclusive

#export OMP_SCHEDULE=STATIC
#export OMP_DISPLAY_ENV=true
#export OMP_NUM_THREADS=48
#export OMP_PROC_BIND=spread
#export OMP_PLACES="{0:1}:48:1"
#export OMP_PLACES='cores(48)'
#mpiexec --map-by ppr:1:socket:pe=24  --report-bindings ./examples/tutorial1
#ddt --connect mpiexec ./project
mpiexec ./project
#mpiexec --map-by socket --bind-to core  --report-bindings ./heterogeneity.exe
#mpiexec --map-by node --bind-to none --report-bindings ./examples/tutorial1
