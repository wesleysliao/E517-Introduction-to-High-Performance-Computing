#include <mpi.h>
#include <stdio.h>
#include <stddef.h>

typedef struct {
  int max_iter;
  double t0;
  double tf;
  double xmax[12];
  double xmin;
} Pars;
 
int main(int argc, char *argv[])
{
  int myid, numprocs, left, right;
  Pars buffer, buffer2;
  MPI_Request request;
  MPI_Status status;
 
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  MPI_Datatype MPI_Pars;
  int count = 5;
  int blocklengths[5] = {1, 1, 1, 12, 1};
  MPI_Aint offsets[5];
    offsets[0] = offsetof(Pars, max_iter);
    offsets[1] = offsetof(Pars, t0);
    offsets[2] = offsetof(Pars, tf);
    offsets[3] = offsetof(Pars, xmax);
    offsets[4] = offsetof(Pars, xmin);
  MPI_Datatype types[5] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

  MPI_Type_create_struct(
    count,
    blocklengths,
    offsets,
    types,    
    &MPI_Pars);
  MPI_Type_commit(&MPI_Pars);
 
  right = (myid + 1) % numprocs;
  left = myid - 1;
  if (left < 0)
    left = numprocs - 1;

  // send myid to the left
  buffer.max_iter = myid;
 
  MPI_Sendrecv(&buffer, 1, MPI_Pars, left, 123, 
     &buffer2, 1, MPI_Pars, right, 123, MPI_COMM_WORLD, &status);

  printf(" Process %d received %d\n",myid,buffer2.max_iter);

  MPI_Finalize();
  return 0;
}


