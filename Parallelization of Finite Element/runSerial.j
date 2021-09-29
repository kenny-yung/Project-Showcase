#!/usr/bin/env bash
# Lecture queue
#SBATCH --account=lect0051
#Outputs of the job
#SBATCH --output=SerialOutput.%j
#SBATCH --error=SerialError.%j
# Wall clock limit
#SBATCH --time=1:00:00

# run the process
../<path to your serial build folder>/2d_Unsteady ./settings.<coarse,medium,fine>.in
