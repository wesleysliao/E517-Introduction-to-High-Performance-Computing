#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

int main(int argc,char **argv) {
  MPI_Init(&argc,&argv);
  int rank,p,i, root = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&p);

  // Make the local vector size constant
  int global_vector_size = 10000;

  double pi = 4.0*atan(1.0);
  
  // initialize the vectors
  double *a, *b;
  a = (double *) malloc(
       global_vector_size*sizeof(double));
  b = (double *) malloc(
       global_vector_size*sizeof(double));
  for (i=0;i<global_vector_size;i++) {
    a[i] = sqrt(i); 
    b[i] = sqrt(i); 
  }

  // compute the dot product
  double sum = 0.0;
  for (i=0;i<global_vector_size;i++) {
    sum += a[i]*b[i];
  }

  if ( rank == root ) {
    printf("The dot product is %g.  Answer should be: %g\n",
           sum,0.5*global_vector_size*(global_vector_size-1));
  }

  free(a);
  free(b);
  MPI_Finalize();
  return 0;
}


