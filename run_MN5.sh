#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --qos=gp_resa
#SBATCH -t 12:00:00
#SBATCH --account=cns119
#SBATCH -o output-%j
#SBATCH -e error-%j
#SBATCH --exclusive


export OMP_DISPLAY_ENV=true
#export OMP_NUM_THREADS=48
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

module purge
module load gcc/13.2.0 openmpi/4.1.5-gcc ddt



#ddt --connect mpirun -n 2  ./heterogeneity.exe ./config/og.xml
srun --nodes=4 --ntasks-per-node=1 --cpus-per-task=112 ./heterogeneity.exe ./config/real_tumour.xml
#ddt --connect mpirun -n 2 ./heterogeneity.exe
#ddt --connect mpirun -n 4 ./heterogeneity.exe
#srun ./capVoxels 2000 1
#srun --nodes=4 --ntasks-per-node=1 --cpus-per-task=112 ./heterogeneity.exe
#./dirichlet_test $voxels 1> ./dirichlet/${voxels}_112_th.log 2>  ./dirichlet/${voxels}_112_th.log

