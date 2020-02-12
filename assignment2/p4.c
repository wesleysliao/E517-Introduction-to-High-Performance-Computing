#include <stdio.h>
#include <stdlib.h>
// link with -lm at compile time
#include <math.h>
#include <omp.h>


int main(int argc, char* argv[]) {
  
  int threads = atoi(argv[1]);
  int const N=1000;
  int i, j, ij;

  double A[N*N];
  double x[N], b[N];

  // initialize the matrix and the vector
  #pragma omp parallel for num_threads(threads) 
  for (ij = 0; ij < N*N; ij++) {
    A[ij] = sin(0.01*(ij));
  }
  
  #pragma omp parallel for num_threads(threads) 
  for (i = 0; i < N; i++) {
    b[i] = cos(0.01*i);
    x[i] = 0.0;
  }

  // matrix vector multiplication
  #pragma omp parallel for num_threads(threads) 
  for (ij = 0; ij < N*N; ij++) {
    #pragma omp atomic
    x[(ij - (ij%N))/N] = x[(ij - (ij%N))/N] + (A[ij]*b[ij%N]);
  }


  printf("x[%d] = %g\n", 505, x[505]);

  return 0;
}
