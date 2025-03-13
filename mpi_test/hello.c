#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv) {
   int rank, size, ierr, plen;
   char pname[MPI_MAX_PROCESSOR_NAME];

   ierr = MPI_Init(&argc, &argv);
   ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   ierr = MPI_Get_processor_name(pname, &plen);

   printf("Hello World, from rank %d of %d procs, host %s\n", rank, size,
          pname);
   ierr = MPI_Finalize();
   return 0;
}
