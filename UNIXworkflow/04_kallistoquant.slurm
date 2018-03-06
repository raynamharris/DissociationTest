#!/bin/bash
#----------------------------------------------------
# Example Slurm job script
# for TACC Stampede2 KNL nodes
#
#   *** Hybrid Job on Normal Queue ***
# 
#       This sample script specifies:
#         10 nodes (capital N)
#         40 total MPI tasks (lower case n); this is 4 tasks/node
#         16 OpenMP threads per MPI task (64 threads per node)
#
# Last revised: 20 Oct 2017
#
# Notes:
#
#   -- Launch this script by executing
#      "sbatch knl.hybrid.slurm" on Stampede2 login node.
#
#   -- Use ibrun to launch MPI codes on TACC systems.
#      Do not use mpirun or mpiexec.
#
#   -- In most cases it's best to specify no more 
#      than 64-68 MPI ranks or independent processes 
#      per node, and 1-2 threads/core. 
#
#   -- If you're running out of memory, try running
#      fewer tasks and/or threads per node to give each 
#      process access to more memory.
#
#   -- IMPI and MVAPICH2 both do sensible process pinning by default.
#
#----------------------------------------------------

#SBATCH -J 04_kallistoquant           # Job name
#SBATCH -o 04_kallistoquant.o%j       # Name of stdout output file
#SBATCH -e 04_kallistoquant.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 4              # Total # of nodes 
#SBATCH -n 14              # Total # of mpi tasks
#SBATCH -t 01:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=rayna.harris@gmail.com
#SBATCH --mail-type=END    # Send email at begin and end of job
#SBATCH -A NeuroEthoEvoDevo       # Allocation name (req'd if you have more than 1)

echo " "
echo " starting 04_kallistoquant Job"
echo " "


cd /scratch/02189/rmharris/DissociationTest/02_filtrimmedreads

module load hdf5/1.8.16
module load zlib/1.2.8
module load kallisto

ibrun ./04_kallistoquant.exe

echo " "
echo " 04_kallistoquant Job Complete"
echo " "
