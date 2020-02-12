#include <stdlib.h> //srand()
#include <time.h>   //time()

#include "mpi.h"

#include "allones.h"

time_t starttime, endtime;

void main() {
  time(&starttime);
  srand(starttime);
  printf("RAND_MAX: %d seed: %d\n", RAND_MAX, starttime);

  time(&endtime);
  printf("execution time: %f s\n", difftime(endtime, starttime));
}
