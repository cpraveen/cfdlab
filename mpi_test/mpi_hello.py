'''
Run like this
    mpirun -np 4 python mpi_hello.py
'''
from mpi4py import MPI
comm = MPI.COMM_WORLD
host = MPI.Get_processor_name()
print("Hello World from host =",host,", rank =",comm.rank)
