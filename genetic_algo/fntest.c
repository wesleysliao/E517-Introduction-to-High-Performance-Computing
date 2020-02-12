#include <stdint.h> //uint8_t
#include <stdio.h>  //sizeof()
#include <stdlib.h> //srand()
#include <time.h>   //time()

#include "mpi.h"

#include "allones.h"

uint8_t testrange[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
double testprobs[10] = {.1, .1, .1, .1, .1, .1, .1, .1, .1, .1};

struct individual indvA = {.genome = {1,1,1,1,1,1,1,1,1,1},
                           .test_scores = {0.0, 0.0, 0.0, 0.0},
                           .tests_done = 0,
                           .final_score = 0.0};

struct individual indvB = {.genome = {2,2,2,2,2,2,2,2,2,2},
                           .test_scores = {0.0, 0.0, 0.0, 0.0},
                           .tests_done = 0,
                           .final_score = 0.0};


time_t starttime, endtime;

void main() {
  time(&starttime);
  srand(starttime);
  printf("RAND_MAX: %d seed: %d\n", RAND_MAX, starttime);

  for(int i = 0; i < 10; i++) {
    printf("%d ", testrange[i]);
  }
  printf("\n");
 
  // Discrete pdf choice
  for(int i = 0; i < 10; i++) {
    printf("%d ", testrange[discrete_pdf_choice(testprobs, 10)]);
  }
  printf("\n");
 
  // array shuffle
  fisher_yates_shuffle(testrange, sizeof(uint8_t), 10);

  for(int i = 0; i < 10; i++) {
    printf("%d ", testrange[i]);
  }
  printf("\n");
  printf("\n");
  
  //2d matrix allocation
  int rows = 5;
  int cols = 4;
  uint8_t** matrix = (uint8_t**) alloc_2d(rows, cols, sizeof(uint8_t));

  for(int i = 0; i < rows; i ++) {
    for(int j = 0; j < cols; j++) {
      printf("%d ", matrix[i][j]);
    }
    printf("\n");
  }

  matrix[1][2] = 1;
  matrix[2][0] = 2;
  matrix[3][1] = 3;
  matrix[0][3] = 4;

  printf("\n");
  for(int i = 0; i < rows; i ++) {
    for(int j = 0; j < cols; j++) {
      printf("%d ", matrix[i][j]);
    }
    printf("\n");
  }

  free(matrix);

  // testgroup randomization
  printf("\n");

  int pop = 10;
  int gs = 2;
  int num_groups = pop/gs;
  int** groups = testgroups(pop, gs);

  printf("\n");
  for(int i = 0; i < num_groups; i ++) {
    for(int j = 0; j < gs; j++) {
      printf("%d ", groups[i][j]);
    }
    printf("\n");
  }
  free(groups);

  //child generation
  struct individual* indvs = malloc(2*sizeof(struct individual));
  indvs[0] = indvA;
  indvs[1] = indvB;
  struct individual child = reproduce(indvs, 2, GENOME_LEN);

  printf("\n");
  for(int i = 0; i < 10; i++) {
    printf("%d ", child.genome[i]);
  }
  printf("\n");

  //genome randomization
  randomize_genomes(indvs, 2, GENOME_LEN);

  printf("\n");
  for(int i = 0; i < 10; i++) {
    printf("%d ", indvs[0].genome[i]);
  }
  printf("\n");


  printf("\n");
  for(int i = 0; i < 10; i++) {
    printf("%d ", indvs[1].genome[i]);
  }
  printf("\n");

  //mutation
  mutate(&indvA, GENOME_LEN, 0.05); 
  printf("\n");
  for(int i = 0; i < 10; i++) {
    printf("%d ", indvA.genome[i]);
  }
  printf("\n");
 

  //test
  int group[2] = {0, 1};
  printf("%d\n", indvs[0].tests_done);
  test(group, 2, indvs, 2, GENOME_LEN);
  printf("%d\n", indvs[0].tests_done);
  test(group, 2, indvs, 2, GENOME_LEN);
  printf("%d\n", indvs[0].tests_done);

  printf("%f\n", indvs[0].test_scores[0]);
  printf("%f\n", indvs[1].test_scores[0]);
  
  //final score
  combine_test_scores(&indvs[0]);
  combine_test_scores(&indvs[1]);

  printf("%f\n", indvs[0].final_score);
  printf("%f\n", indvs[1].final_score);
  
  printf("\n");
  
  //convergence
  struct individual kiddos[2];
  kiddos[0] = reproduce(indvs, 2, 10);
  kiddos[1] = reproduce(indvs, 2, 10);

  kiddos[0] = indvs[0];
  printf("converged: %d\n", population_converged(indvs, kiddos, 2, 10));
  printf("converged: %d\n", population_converged_diff(indvs, kiddos, 2, 10, 0.1));
  indvs[1] = indvs[0];
  kiddos[0] = indvs[0];
  kiddos[1] = indvs[0];
  printf("converged: %d\n", population_converged(indvs, kiddos, 2, 10));
  printf("converged: %d\n", population_converged_diff(indvs, kiddos, 2, 10, 0.1));

  printf("\n");
  time(&endtime);
  printf("execution time: %fs\n", difftime(endtime, starttime));

  struct run_conditions paramset = {.runs_to_aggr = 10 ,
                                    .generations_max = 20,
                                    .genome_size = 10,
                                    .population_size = 10,
                                    .testgroup_size = 1,
                                    .tests_per_indv = 1,
                                    .parents_per_repro = 2,
                                    .children_per_repro = 1,
                                    .mutation_rate = 0.01};

  struct run_data single_run = genesis_run(GENERATIONS_MAX,
                                           POPULATION_SIZE,
                                           GENOME_LEN, 1, 1, 2, 1, 0.01);

  printf("converged: %d\n", single_run.converged);
  
}
