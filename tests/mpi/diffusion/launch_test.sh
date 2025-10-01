#!/bin/bash
#SBATCH --nodes=50
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --qos=gp_resa
#SBATCH -t 2:00:00
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

srun --nodes=50 --ntasks-per-node=1 --cpus-per-task=112 ./diffusion ./config/100nodes.xml >> large.log
srun --nodes=50 --ntasks-per-node=1 --cpus-per-task=112 ./diffusion >> small.log
#ddt --connect mpirun -n 100 ./diffusion ./config/100nodes.xml 
#ddt --connect mpirun -n 2 ./cell_snd
