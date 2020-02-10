#!/bin/sh
#PBS -V
#PBS -N mpi_test
#PBS -o pbs_out.txt
#PBS -e pbs_error.txt
#PBS -l nodes=1:ppn=4
#PBS -q test
cd $PBS_O_WORKDIR
echo "Job began"
date
mpiexec ./hello
echo "Job ended"
date
