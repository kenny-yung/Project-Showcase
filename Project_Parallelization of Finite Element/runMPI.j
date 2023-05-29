#!/usr/bin/env bash
# Lecture queue
#SBATCH --account=lect0051
### Select the number of MPI PEs to use. 
#SBATCH --ntasks=<1,2,4,6,8,12,16>
#Outputs of the job
#SBATCH --output=MPIOutput.%j
#SBATCH --error=MPIError.%j
# Wall clock limit
#SBATCH --time=1:00:00
### run the process
### select a mesh
$MPIEXEC $FLAGS_MPI_BATCH ../<path to your serial build folder>/2d_Unsteady ./settings.<coarse,medium,fine>.in
