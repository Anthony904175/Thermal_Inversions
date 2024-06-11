#!/bin/bash
#SBATCH --mail-user=aosborne3@ucmerced.edu  
#SBATCH --mail-type=ALL  #uncomment the first two lines if you want to receive the email notifications
#SBATCH --job-name=octoPIB  # name of the job
#SBATCH --partition=pi.apribram-jones  # partition to be used (defq, gpu or intel)
#SBATCH --time=01:00:00    # walltime (up to 96 hours)
#SBATCH --nodes=1      # number of nodes
#SBATCH --ntasks=1     # number of tasks (i.e. parallel processes) to be started
#SBATCH --cpus-per-task=56
#SBATCH --output=pib.stdout
#SBATCH --exclusive
#SBATCH --export=ALL

module purge
module load mvapich2/2.3.6-intel-2021.4.0
module load mpi/2021.4.0
module load mkl/2021.4.0
module load gsl
module load gcc/11.2.0

export TEMPDIRPATH=$SCRATCH/tmp
OCTOPUS=/home/aosborne3/user-modules/octopus-12.1/src/octopus 

export OMP_NUM_THREADS=1
mpirun -n 56 $OCTOPUS &> pib.out
