#include <stdint.h> //uint8_t
#include <stdlib.h> //rand

#define GENERATIONS_MAX         5000
#define POPULATION_SIZE         10

#define GENOME_LEN              10
#define GENOME_TYPE             uint8_t
#define TESTS_PER_INDV          4
#define PARENTS_PER_REPRO       2


#include "gatypes.h"

GENOME_TYPE random_gene() {
  return (rand() % 2);
}

GENOME_TYPE mutate_gene(GENOME_TYPE original_gene) {
  return !original_gene;
}

void test(int* group, int testgroup_size, struct individual* individuals, int population_size, int genome_size) {
  
  for (int member_ndx = 0; member_ndx < testgroup_size; member_ndx++) {
    double gene_total = 0;
    for (int gene_ndx = 0; gene_ndx < genome_size; gene_ndx++) {
      gene_total += individuals[group[member_ndx]].genome[gene_ndx]; 
    }
    individuals[group[member_ndx]].test_scores[individuals[group[member_ndx]].tests_done++] = gene_total;
  }
}

void combine_test_scores(struct individual* indv) {
  indv->final_score = 0.0;
  for (int test_ndx = 0; test_ndx < indv->tests_done; test_ndx++) {
    indv->final_score += indv->test_scores[test_ndx];
  }
}

#include "ga.h"
