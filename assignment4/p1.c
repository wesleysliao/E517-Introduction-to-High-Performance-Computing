#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
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
 
  time_t starttime;
  if (myid == 0)
    time(&starttime);
  MPI_Barrier(MPI_COMM_WORLD);

 
  long divisions = strtol(argv[1], NULL, 10);
  const double radius = 2.0;
  const double div_width = radius/divisions;

  double mysum = 0.0;

  for(long division = myid; division < divisions; division += numprocs) {
    double l = (division+0.5)*div_width;
    double h = sqrt((radius*radius)-(l*l));
    
    mysum += h * div_width;
  }

  double total = 0.0;
  MPI_Reduce( &mysum, &total, 1,
              MPI_DOUBLE, MPI_SUM,
              0, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);

  time_t endtime;
  if(myid == 0){
    time(&endtime);
    double timediff = difftime(starttime, endtime);

    printf("the total is %1.15f\n", total);
    printf("true value:  %1.15f\n", M_PI);
    printf("req. accuracy:           ^\n");
    printf("%20d, %d\n", starttime, endtime);
    printf("execution time (wall clock): %1.20f s\n", timediff);
  }

  MPI_Finalize();
  return 0;
}
