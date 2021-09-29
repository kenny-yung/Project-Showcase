#!/usr/bin/env bash
# Lecture queue
#SBATCH --account=lect0051
#Outputs of the job
#SBATCH --output=OMPOutput.%j
#SBATCH --error=OMPError.%j
# Wall clock limit
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
## Specify the desired number of threads
#SBATCH --cpus-per-task=<1,2,4,6,8,12>
# run the process
../<path to your OMP build folder>/2d_Unsteady ./settings.<coarse,medium,fine>.in
