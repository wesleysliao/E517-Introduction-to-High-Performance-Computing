#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

int main(int argc,char **argv) {
  MPI_Init(&argc,&argv);
  int rank,p,i, root = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&p);

  // Make the global vector size constant
  int global_vector_size = 10000;

  int local_vector_size = (global_vector_size / p);

  double pi = 4.0*atan(1.0);
  
  // initialize the vectors
  double *a, *b;
  a = (double *) malloc(
       local_vector_size*sizeof(double));
  b = (double *) malloc(
       local_vector_size*sizeof(double));
  for (i=0;i<local_vector_size;i++) {
    a[i] = sqrt(i+(rank*local_vector_size)); 
    b[i] = sqrt(i+(rank*local_vector_size)); 
  }

  double mysum = 0.0;
  for(i = 0; i < local_vector_size; i++)
  {
    mysum += a[i]*b[i];
  }

  // compute the dot product
  double total = 0.0;
  MPI_Reduce(
    &mysum,
    &total,
    1,
    MPI_DOUBLE,
    MPI_SUM,
    0,
    MPI_COMM_WORLD);

  if ( rank == 0) {
    printf("The dot product is %g.  Answer should be: %g\n",
           total, 0.5*global_vector_size*(global_vector_size-1));
  }

  free(a);
  free(b);
  MPI_Finalize();
  return 0;
}


