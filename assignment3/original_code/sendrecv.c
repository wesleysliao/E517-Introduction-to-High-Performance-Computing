#include "mpi.h"
#include <stdio.h>
 
int main(int argc, char *argv[])
{
    int myid, numprocs, left, right;
    int buffer, buffer2;
    MPI_Request request;
    MPI_Status status;
 
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
 
    right = (myid + 1) % numprocs;
    left = myid - 1;
    if (left < 0)
        left = numprocs - 1;

    // send myid to the left
    buffer = myid;
 
    MPI_Sendrecv(&buffer, 1, MPI_INT, left, 123, 
       &buffer2, 1, MPI_INT, right, 123, MPI_COMM_WORLD, &status);

    printf(" Process %d received %d\n",myid,buffer2);

    MPI_Finalize();
    return 0;
}


