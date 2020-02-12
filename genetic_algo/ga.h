#ifndef GENETIC_ALGO_ENGINE_H
#define GENETIC_ALGO_ENGINE_H

#include <math.h> //fabs()
#include <stdint.h> //uint8_t
#include <stdio.h> //file io
#include <stdlib.h> //rand()
#include <string.h> //memcpy()


int compare_int(const void* int_a, const void* int_b) {
  if(*(int*)int_a < *(int*)int_b) return -1;
  if(*(int*)int_a == *(int*)int_b) return 0;
  if(*(int*)int_a > *(int*)int_b) return 1;
}

//
// Memory and array functions
//

void memswap(void* dest, void* src, size_t len) {
  void* buffer = malloc(len);

  memcpy(buffer, dest, len);
  memcpy(dest, src, len);
  memcpy(src, buffer, len);

  free(buffer);
}
 
void** alloc_2d(int num_rows, int num_cols, size_t typesize) {
  void* values = calloc(num_rows * num_cols, typesize); 
  void** rows = malloc(num_rows * sizeof(void*));
  for(int i = 0; i < num_rows; i++) {
    rows[i] = values + (i * num_cols * typesize);
  }

  return rows;
}

//
// Random functions
//

int uniform_choice(int len){
  return (rand() % len);
}

int discrete_pdf_choice(double *probabilities, int len) {
  // Selects one index from array according to given probabilities
  // probabilities should sum to 1.

  int rand_value = rand();
  double cumulative_prob = 0.0;

  for(int i = 0; i < (len - 1); i++) {
    cumulative_prob += probabilities[i];
    if ( rand_value <= (int) (cumulative_prob * RAND_MAX)) {
      return i;
    }
  }
  return (len - 1);
}

void fisher_yates_shuffle(void *array, size_t typesize, unsigned int len) {
  // Shuffles elements of an array
  for(int i = 0; i < (len - 2); i++) {
    int j = (rand() % (len - i)) + i;
    
    memswap(array + (typesize * i),array + (typesize * j), typesize);
  }
}

//
// Genetic Algorithm functions
//

void randomize_genomes(struct individual* individuals, 
                       int population_size, 
                       int genome_len) {
  for(int indv_ndx = 0; indv_ndx < population_size; indv_ndx++) {
    for(int gene_ndx = 0; gene_ndx < genome_len; gene_ndx++) {
      individuals[indv_ndx].genome[gene_ndx] = random_gene();
    }
  }
}

int** testgroups(int population_size, int testgroup_size) {
  int* permutation = malloc(population_size * sizeof(int));
  for (int i = 0; i < population_size; i++) {
    permutation[i] = i;
  }
  fisher_yates_shuffle(permutation, sizeof(int), population_size);
    
  int num_groups = population_size / testgroup_size;
  int** groups = malloc(num_groups * sizeof(int*));
  for (int group_ndx = 0; group_ndx < num_groups; group_ndx++) {
    groups[group_ndx] = permutation + (group_ndx * testgroup_size);
  }
  
  return groups;
}


struct individual reproduce(struct individual* parents, 
                            unsigned int num_parents,
                            unsigned int genome_size){

  // Create a range of all gene indices
  int* section_range = malloc((genome_size - 1) * sizeof(int));
  for(int gene_ndx = 1; gene_ndx < genome_size; gene_ndx++) {
    section_range[gene_ndx-1] = gene_ndx;
  }
  //shuffle the gene indices
  fisher_yates_shuffle(section_range, sizeof(int), genome_size - 1);

  //pick the first number of indices as we want sections
  int sections  = 2;
  int* section_start = malloc(sections * sizeof(int));
  section_start[0] = 0;
  for(int section_ndx = 1; section_ndx < sections; section_ndx++) {
    section_start[section_ndx] = section_range[section_ndx-1];
  }
  qsort(section_start, sections, sizeof(int), compare_int);

  struct individual child;
  child.tests_done = 0;
  child.final_score = 0.0;
  memcpy(child.parents, parents, (sizeof(unsigned int)*num_parents));
  
  //copy each section to child, taking turns from parents
  for(int section_ndx = 0; section_ndx < sections; section_ndx++) {
    int section_size;
    if(section_ndx == sections - 1) {
      section_size = genome_size - section_start[section_ndx];
    } else {
      section_size = section_start[section_ndx + 1] - section_start[section_ndx];
    }
     
    int parent_ndx = section_ndx % num_parents;
    memcpy(&child.genome[section_start[section_ndx]], 
           &parents[parent_ndx].genome[section_start[section_ndx]],
           section_size * sizeof(GENOME_TYPE));
  }
  return child;
}

void mutate(struct individual* indv, int genome_size, double mutation_rate) {
  int mutate_index = rand() % (int)(1/mutation_rate);

  if(mutate_index < genome_size) {
    indv->genome[mutate_index] = mutate_gene(indv->genome[mutate_index]);
  }
}

int population_converged(struct individual* parents, 
                         struct individual* children,
                         int population_size,
                         int genome_size){
  
  for(int gene_ndx = 0; gene_ndx < genome_size; gene_ndx++) {
    for(int indv_ndx = 0; indv_ndx < population_size; indv_ndx++) {
      if(parents[indv_ndx].genome[gene_ndx]
         != children[indv_ndx].genome[gene_ndx]) {
        return 0;
      }
      for(int comp_ndx = 0; comp_ndx < population_size; comp_ndx++) {
        if( comp_ndx == indv_ndx ) {
          continue;
        }
        if(children[indv_ndx].genome[gene_ndx]
           != children[comp_ndx].genome[gene_ndx]) {
          return 0;
        }
      }
    }
  }
  return 1;
}

int population_converged_diff(struct individual* parents, 
                              struct individual* children,
                              int population_size,
                              int genome_size,
                              double epsilon){
  
  for (int gene_ndx = 0; gene_ndx < genome_size; gene_ndx++) {
    for (int indv_ndx = 0; indv_ndx < population_size; indv_ndx++) {
      if (fabs(parents[indv_ndx].genome[gene_ndx]
               - children[indv_ndx].genome[gene_ndx])
          >= epsilon) {
        return 0;
      }
      for (int comp_ndx = 0; comp_ndx < population_size; comp_ndx++) {
        if ( comp_ndx == indv_ndx ) continue;
       
        if (fabs(children[indv_ndx].genome[gene_ndx]
                 - children[comp_ndx].genome[gene_ndx])
            >= epsilon) {
          return 0;
        }
      }
    }
  }
  return 1;
}

//
// Main algorithm
//

struct run_data genesis_run(int generations_max,
                 int population_size,
                 int genome_size,
                 int testgroup_size,
                 int tests_per_indv,
                 int parents_per_repro,
                 int children_per_repro,
                 double mutation_rate) {

  struct run_data data = {0};
  randomize_genomes(data.generations[0].individuals, population_size, genome_size); 

  for(int generation = 0; generation < generations_max; generation++) {
    printf("%d \n", generation);


    struct individual* individuals = data.generations[generation].individuals;

    //////
    // Test each group
    for(int test_ndx = 0; test_ndx < tests_per_indv; test_ndx++) {

      int** groups = testgroups(population_size, testgroup_size);
      int num_groups = population_size / testgroup_size;

      for(int group_ndx = 0; group_ndx < num_groups; group_ndx++) {
        test(groups[group_ndx],
             testgroup_size,
             individuals,
             population_size,
             genome_size);
      }
    }
    //////


    //////
    // Combine scores and calculate mating proability
    double* mating_probability = malloc(population_size*sizeof(double));
    double prob_sum = 0.0;
    for(int indv_ndx = 0; indv_ndx < population_size; indv_ndx++) {
      combine_test_scores(&individuals[indv_ndx]); 
      
      mating_probability[indv_ndx] = individuals[indv_ndx].final_score; 
      prob_sum += mating_probability[indv_ndx];
    }

    // normalize
    for(int prob_ndx = 0; prob_ndx < population_size; prob_ndx++) {
      mating_probability[prob_ndx] /= prob_sum;
    }
    //////

    
    //////
    // Reproduce
    if(generation < (generations_max - 1)){ 
      
      struct individual* children = data.generations[generation+1].individuals;
      
      for (int kid_ndx = 0; 
           kid_ndx < population_size; 
           kid_ndx += children_per_repro) {

        struct individual** parents = malloc(parents_per_repro
                                             * sizeof(struct individual*));
        int* parent_indices = malloc(parents_per_repro*sizeof(int));
        for (int parent_ndx = 0; 
             parent_ndx < parents_per_repro;
             parent_ndx++) {

          int same_parent;
          do {
            //select a parent
            parent_indices[parent_ndx] = discrete_pdf_choice(mating_probability,
                                                             population_size);
            parents[parent_ndx] = &individuals[parent_indices[parent_ndx]];
            
            //ensure they haven't already been selected
            same_parent = 0;
            for (int check_ndx = 0; check_ndx < parent_ndx; check_ndx++) {
              same_parent += (parent_indices[parent_ndx] 
                              == parent_indices[check_ndx]);
            }
          } while(same_parent > 0);
        }

        children[kid_ndx] = reproduce(*parents,
                                      parents_per_repro,
                                      genome_size);

        mutate(&children[kid_ndx], genome_size, mutation_rate);

        for(int i = 0; i < 10; i++) printf("%d ", children[kid_ndx].genome[i]); printf("\n");

        free(parents);
        free(parent_indices);
      }


      if(population_converged(individuals, 
                              children,
                              population_size, 
                              genome_size)) {
        data.converged = 1;
        //break; 
      }
    }
    free(mating_probability);
    //////

  }
}


void aggregated_runs(struct run_conditions params, FILE* outfile) {

  fwrite(&params, sizeof(struct run_conditions), 1, outfile);

  struct run_data* run_record = malloc(sizeof(struct run_data));
  
  for(int run = 0; run < params.runs_to_aggr; run++){
    run_record[run] = genesis_run(params.generations_max,
                                  params.population_size,
                                  params.genome_size,
                                  params.testgroup_size,
                                  params.tests_per_indv,
                                  params.parents_per_repro,
                                  params.children_per_repro,
                                  params.mutation_rate);

    fwrite(&run_record[run], sizeof(struct run_data), 1, outfile);
  }
  
  free(run_record);
}

#endif
