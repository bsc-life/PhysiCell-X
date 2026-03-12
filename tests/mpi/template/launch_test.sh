#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --qos=gp_debug
#SBATCH -t 02:00:00
#SBATCH --account=cns119
#SBATCH -o output-%j
#SBATCH -e error-%j
#SBATCH --exclusive


export OMP_DISPLAY_ENV=false
#export OMP_NUM_THREADS=48
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

module purge
module load gcc/13.2.0 openmpi/4.1.5-gcc ddt

#export UCX_TLS=rc,t

#srun ./cell_snd 100
#ddt --connect mpirun -n 2 ./cell_snd 100
ddt --connect mpirun -n 2 ./cell_snd
