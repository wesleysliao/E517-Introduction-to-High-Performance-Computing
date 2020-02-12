// Copyright 2019 Wesley Liao
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include <string.h>

#include "mmio.h"


int main(int argc, char *argv[]) {

  //
  // Open the matrix file
  //
  char *matrix_filename = argv[1];

  FILE *matrix_file;
  matrix_file = fopen(matrix_filename,"r");

  if (matrix_file == 0) {

    printf(strcat(strcat("could not open file \"", matrix_filename), "\"\n"));

  }
  else {

    //
    // Parse matrix data from file
    //
      
    MM_typecode matrix_typecode;
    int number_of_rows;
    int number_of_columns;
    int number_of_nonzeros;
      
    mm_read_banner(matrix_file, &matrix_typecode);
    printf(strcat(mm_typecode_to_str(matrix_typecode), "\n"));

    if (mm_is_coordinate(matrix_typecode) 
        && mm_is_pattern(matrix_typecode) 
        && mm_is_symmetric(matrix_typecode)) {
        
      mm_read_mtx_crd_size(matrix_file, &number_of_rows, &number_of_columns, &number_of_nonzeros);

      printf("%i %i %i\n", number_of_rows, number_of_columns, number_of_nonzeros);

      int *matrix_data = calloc(number_of_rows*number_of_columns, 4);

      for (int i = 0; i<number_of_nonzeros; i++){
        // Read in the next coordinate pair from the matrix file
        int row;
        int col;
        fscanf(matrix_file, "%i %i", &row, &col);
        printf("%i %i\n", row, col);

        // Decrement the indexes because MM format indexes from 1 instead of 0
        row -= 1;
        col -= 1;

        // Set the corresponding matrix location to 1
        *(matrix_data + ((row*number_of_rows) + col)) = 1;
        // Because the matrix is symmetric, also set the opposite location
        *(matrix_data + ((col*number_of_rows) + row)) = 1;
      }

      // Print our imported matrix
      for (int i = 0; i < (number_of_rows*number_of_columns); i++) {
        // Because our data is 1d, print one character at a time
        // and when we hit a multiple of the width print a newline
        if (i % number_of_rows == 0)
            printf("\n");
        printf("%i", matrix_data[i]);
      }
      printf("\n");

      //
      // Create our vector
      //

      double *operation_vector = calloc(number_of_columns, 8);

      for (int i = 0; i < number_of_columns; i++) {
        *(operation_vector+i) = sin((2*M_PI*i) / (number_of_columns-1));
      }

      //
      // Compute matrix-vector product
      //

      double *result_vector = calloc(number_of_rows, 8);

      for (int row = 0; row < number_of_rows; row++) {
        for (int col = 0; col < number_of_columns; col++) {

            int matrix_value = *(matrix_data + ((row*number_of_rows) + col));

            if ( matrix_value == 1 ) { 
              *(result_vector + col) += *(operation_vector + col);
            }
        }
      }

      for (int i = 0; i < number_of_columns; i++) {
        printf("%f, %f\n", *(operation_vector+i), *(result_vector+i));
      }
      
      //
      // Calculate l1 norm of result vector
      //

      double l1norm = 0.0;
      for (int i = 0; i < number_of_rows; i++) {
          l1norm += *(result_vector + i);
      }
      printf("l1 norm of result: %f\n", l1norm);

      //
      // Calculate total number of floating-point operations
      //

      int flops = 
        (3*number_of_columns)// 3 for each element in our operation vector
        + number_of_nonzeros // one fp addition per non-zero value in the matrix
        + number_of_columns; // final l1 norm sum

      printf("total number of flops: %i\n", flops);

      // Clean up allocated memory
      free(matrix_data);
      free(operation_vector);
      free(result_vector);
    }
    else {
        printf("bad matrix type. abort.\n");
    }
     
    fclose(matrix_file);
  }
}
