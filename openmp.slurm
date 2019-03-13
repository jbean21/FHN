#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3882mb
#SBATCH --time=00:01:00
#PBS -q px425

# Change into the directory where I submitted the job
cd $SLURM_SUBMIT_DIR 

# Load the modules used when compiling my code
module load intel impi

# Set the number of threads to use, this should match
# the number of processors per node (ppn) above.
export OMP_NUM_THREADS=16

# Execute gravity executable
./a.out -N 1

exit 0
