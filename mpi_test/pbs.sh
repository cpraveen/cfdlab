#!/bin/sh
#PBS -V
#PBS -m ae
#PBS -M cpraveen@gmail.com
#PBS -N mpi_test
#PBS -o pbs_out.txt
#PBS -e pbs_error.txt
#PBS -l nodes=1:ppn=4
#PBS -q test
cd $PBS_O_WORKDIR
echo "Job began"
echo "Working directory is $PBS_O_WORKDIR"

# Calculate the number of processors allocated to this run.
NPROCS=`wc -l < $PBS_NODEFILE`

# Calculate the number of nodes allocated.
NNODES=`uniq $PBS_NODEFILE | wc -l`

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors across ${NNODES} nodes

# OpenMPI will automatically launch processes on all allocated nodes.
# MPIRUN=`which mpirun`
# ${MPIRUN} -machinefile $PBS_NODEFILE -np ${NPROCS} my-openmpi-program

mpiexec ./hello

echo "Job ended"
echo Time is `date`
