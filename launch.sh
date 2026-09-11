#!/bin/bash -l
# The -l above is required to get the full environment with modules

#SBATCH --job-name=interp
#SBATCH --partition=compute
#SBATCH --time=0-02:00:00
#SBATCH --ntasks=256
#SBATCH --mem-per-cpu=3500M
#SBATCH --mail-type=BEGIN,FAIL,END

srun --mpi=pmix ./interpolate_MPI.x > log.out 2>err.out
