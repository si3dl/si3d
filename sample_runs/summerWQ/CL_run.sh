#!/bin/bash -l
# NOTE the -l flag!

# If you need any help, please email help@cse.ucdavis.edu

# Name of the job 
#SBATCH -J CLWQ

# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o CLWQ-%j.output
#SBATCH -e CLWQ-%j.output
#SBATCH --exclude=agate-8,agate-6,agate-5,agate-7,agate-45,agate-0


# no -n here, the user is expected to provide that on the command line.

# The useful part of your job goes below

# run one thread for each one the user asks the queue for
# hostname is just for debugging
hostname
export OMP_NUM_THREADS=$SLURM_NTASKS
module load benchmarks intel

# The main job executable to run: note the use of srun before it
time srun psi3d
