#include <stdio.h>
#include <omp.h>

// compute the dot product of two vectors

int main() {
  int const N=100;
  int i;
  double a[N], b[N];
  double dot_prod = 0.0;

  int thread_id;
  //Arbitrarily initialize vectors a and b
  for (i = 0; i < N; i++) {
    a[i] = 3.14;
    b[i] = 6.67;
  }

  #pragma omp parallel private(thread_id)
  {
    thread_id = omp_get_thread_num();
    printf("This thread is: %d\n", thread_id);
    #pragma omp for
    for(i=0; i < N; i++) {
      // sum up the element-wise product of the two arrays
      #pragma omp atomic 
      {
        dot_prod = dot_prod + (a[i] * b[i]);
      }
    }
  }

  printf("Dot product of the two vectors is %g\n", dot_prod);
  
  return 0;
}
