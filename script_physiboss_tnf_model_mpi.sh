#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=24
#SBATCH --constraint=hyperthreading
#SBATCH -t 08:00:00
#SBATCH -o output-%j
#SBATCH -e error-%j
#SBATCH --exclusive

#-----------------------------------------------------------------------------------
# Except for OMP_SCHEDULE=STATIC/DYNAMIC, the rest are the same as that for BioFVM_X
# Better to use OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK, this way will only need to
# change --cpus-per-task=<no_of_threads>
#-----------------------------------------------------------------------------------

export OMP_DISPLAY_ENV=true
export OMP_SCHEDULE=STATIC
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

#export OMP_PLACES="{0:1}:48:1"
#export OMP_PLACES='cores(48)'

#--------------------------------------------
# Simplest Execution, can be used for testing
#--------------------------------------------

# mpiexec ./spheroid_TNF_model_mpi
# ddt --connect mpiexec ./spheroid_TNF_model_mpi
# ddt --connect srun ./spheroid_TNF_model_mpi
#-------------------------------------------------------------
# Better to use --map-by ppr syntax when measuring performance
# MN4, best configuration is 1 MPI process per socket and
# 24 OpenMP threads per MPI process.
# (This is because we have 2 sockets with 24 cores each)
#-------------------------------------------------------------

 mpiexec --map-by ppr:1:socket:pe=24  --report-bindings ./spheroid_TNF_model_mpi

#-------------------------------------
# Uncomment if using DDT for debugging (1) if mpiexec doesn't connect then use (2) srun
#-------------------------------------

# ddt --connect mpiexec ./spheroid_TNF_model_mpi
# ddt --connect srun ./spheroid_TNF_model_mpi


#---------------------------------------------
# This is for taking outputs of 100 experiments
#---------------------------------------------

# for i in `seq 1 100`
#         do
#                 mpiexec ./project
#                 mkdir /gpfs/scratch/bsc99/bsc99102/PHYSICELL/PARALLEL_NEW/parallel_exp_$i
#                 ls -d ./output/* > test.txt && xargs -a test.txt cp -t /gpfs/scratch/bsc99/bsc99102/PHYSICELL/PARALLEL_NEW/parallel_exp_$i
#                 rm ./output/*.mat ./output/*.xml ./output/*.svg CELLS_RANK_*
#         done


#----------------------------------------------------------
# Using OpenMPI binding policies (non ppr execution syntax)
#----------------------------------------------------------

#mpiexec --map-by socket 	--bind-to core  --report-bindings ./heterogeneity.exe  <--- This is for heterogeneity example
#mpiexec --map-by node 		--bind-to none 	--report-bindings ./examples/tutorial1 <--- This was for BioFVM
