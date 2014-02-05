#!/usr/bin/env zsh
np=4
#BSUB -J SiScLab	  # job name
##BSUB -U "p_sisclab#2"    # priority queue
#BSUB -o run.qout	  # job output (use %J for job ID)
#BSUB -W 00:30    	  # limits in hours:minutes
#BSUB -M 512              # memory in MB
#BSUB -n 4               # number of slots
#BSUB -a openmpi  	  # MPI environment
#BSUB -m mpi-l            # host group (mpi-s or mpi-l for Westmere) "bcs openmpi"

module load MISC VTK

$MPIEXEC $FLAGS_MPI_BATCH ../src/solver/2d_Unsteady_Diffusion>run_np_superfine_$np.out
