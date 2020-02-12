#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
//#include <stddef.h>
#include <math.h>
#include <time.h>

int main(int argc, char *argv[])
{
  int myid, numprocs;
  MPI_Request request;
  MPI_Status status;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  
  time_t start;
  if(myid == 0){
    time(&start);
  }
  long divisions = strtol(argv[1], NULL, 10);
  const double interval = M_PI/2.0;
  const double div_width = interval/divisions;

  double mysum = 0.0;

  for(long division = myid; division < divisions; division += numprocs) {
    double x = division*div_width;
    double h = cos(x)*sin(2*x);
    
    //printf("rank %d computing div %d of %d: x=%f, h=%f %f\n", myid, division, divisions-1, x, h, div_width);
    mysum += h * div_width;
  }

  double total = 0.0;
  MPI_Reduce( &mysum, &total, 1,
              MPI_DOUBLE, MPI_SUM,
              0, MPI_COMM_WORLD);

  if(myid == 0){
    time_t end;
    time(&end);
    double timediff = difftime(end, start);
    printf("the total is %1.15f\n", total);
    printf("true value   %1.15f\n", 2.0/3.0);
    printf("req. accuracy:       ^\n");
    printf("elapsed time %f s\n", timediff);
  }

  MPI_Finalize();
  return 0;
}
