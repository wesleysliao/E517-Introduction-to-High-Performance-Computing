#ifndef GENETIC_ALGO_ENGINE_TYPES_H
#define GENETIC_ALGO_ENGINE_TYPES_H

//
// Datatypes
//

struct individual {
  GENOME_TYPE genome[GENOME_LEN];
  unsigned int parents[PARENTS_PER_REPRO];

  double test_scores[TESTS_PER_INDV];
  unsigned int tests_done;
  double final_score;
};

struct generation {
  struct individual individuals[POPULATION_SIZE];

  double high_score;
  double mean_score;
  double stdv_score;
};

struct run_conditions {
  unsigned int runs_to_aggr;
  unsigned int generations_max;
  unsigned int genome_size;
  unsigned int population_size;
  unsigned int testgroup_size;
  unsigned int tests_per_indv;
  unsigned int parents_per_repro;
  unsigned int children_per_repro;
  double mutation_rate;
};

struct run_data {
  unsigned int run_count;
  unsigned int generations_completed;
  unsigned int converged;

  struct generation generations[GENERATIONS_MAX];
};

#endif
