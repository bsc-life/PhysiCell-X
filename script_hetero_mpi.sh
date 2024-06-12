#!/bin/bash

#SBATCH --nodes=3
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH -t 02:00:00
#SBATCH -o output-%j
#SBATCH -e error-%j
#SBATCH --exclusive
#SBATCH --qos=debug

# Total MPI processes = 2 nodes x 1 = 2
# Total OpenMP threads = 2 MPI processes x 48 threads = 96
# Total nodes used = 2
# Total cores used = 2 x 48 = 96 = No. of OpenMP threads

export OMP_SCHEDULE=STATIC
export OMP_DISPLAY_ENV=true
export OMP_NUM_THREADS=48		#can be set equal to $SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread

#Simplest execution 

# mpiexec --map-by ppr:1:node:pe=24 ./heterogeneity.exe config/het.xml
mpiexec ./heterogeneity.exe config/het.xml
# Execution using --map-by ppr syntax to create 1 MPI process per node
# and 48 threads per MPI process and also report the MPI+OpenMP bindings in
# the standard error file, use this syntax for performance, possibly best
# is to use 1 MPI process per socket and 24 threads per MPI process
# (also see the TNF MPI example) 

#mpiexec --map-by ppr:1:node:pe=48  --report-bindings ./heterogeneity.exe

# Syntax for connecting to Arm DDT through "srun"

#ddt --connect srun ./heterogeneity.exe


