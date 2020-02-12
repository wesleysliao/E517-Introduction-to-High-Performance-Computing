#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "mpi.h"

#define COLS 1000
#define ROWS 1000
#define TEMP 50.
#define DEBUG 0
#define EPS 1e-6
#define I_FIX COLS/2
#define J_FIX ROWS/2

double max_abs(double** m1, double** m2, int width, int height){
    double max_val = DBL_MIN;
    for (int i = 1; i < height-1; i++)
        for (int j = 1; j < width-1; j++){
            if (fabs(m1[i][j] - m2[i][j]) > max_val) {
                max_val = fabs(m1[i][j] - m2[i][j]);
            }
        }
    return max_val;
}

void print_matrix(double** matrix, int width, int height){
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++)
            printf("%2.0f", matrix[i][j]);
        printf("\n");
    }
}

void print_matrix_csv(double** matrix, int width, int height){
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            printf("%17.15f", matrix[i][j]);
            if(j < width-1)
              printf(",");
            }
        printf("\n");
    }
}

void file_matrix_csv(FILE* file, double** matrix, int width, int height){
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            fprintf(file, "%17.15f", matrix[i][j]);
            if(j < width-1)
              fprintf(file, ",");
            }
        fprintf(file, "\n");
    }
}
void copy_matrix(double** dest, double** source, int width, int height) {
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
            dest[i][j] = source[i][j];
}

void copy_row_to_buf(double** matrix, double* buffer, int width, int row) {
  for (int j = 0; j < width; j++)
    buffer[j] = matrix[row][j+1];
}

void copy_col_to_buf(double** matrix, double* buffer, int height, int col) {
  for (int i = 0; i < height; i++)
    buffer[i] = matrix[i+1][col];
}

void copy_row_from_buf(double** matrix, double* buffer, int width, int row) {
  for (int j = 0; j < width; j++)
     matrix[row][j+1] = buffer[j];
}

void copy_col_from_buf(double** matrix, double* buffer, int height, int col) {
  for (int i = 0; i < height; i++)
    matrix[i+1][col] = buffer[i];
}

double** alloc_matrix(int width, int height){
    double** matrix;
    matrix = (double**) malloc(height * sizeof(double *));
    matrix[0] = (double*) malloc(height * width * sizeof(double));
    for (int i = 1; i < height; i++)
        matrix[i] = matrix[0] + i*width;
    return matrix;
}

void compute_new_values(double** old_matrix, double** new_matrix, int width, int height){
    for (int i = 1; i < height-1; i++)
        for (int j= 1; j < width-1; j++)
            new_matrix[i][j] =
                    0.25 * (old_matrix[i-1][j] + old_matrix[i+1][j] +
                            old_matrix[i][j-1] + old_matrix[i][j+1]);
}

void init_matrix(double** matrix, int width, int height){
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++) {
            matrix[i][j] = 0.;
        }
}

int is_power_of_2(int test_value) {
  if (test_value == 2 || test_value == 1)
    return 1;
  else if (test_value % 2 == 0)
    return is_power_of_2(test_value/2);
  else
    return 0;
}

int root_of_next_square(int value)
{
  int test = 1;
  for(test=1; test*test < value; test++);
  return test;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int mpi_rank, mpi_size, root=0;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    time_t start, end;
    if(mpi_rank == root)
      time(&start);

    
    if(mpi_rank == root && !is_power_of_2(mpi_size)) {
      printf("number of processors must be power of 2. Aborting.\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
      return -1;
    }
    
    int chunk_cols = root_of_next_square(mpi_size);
    int chunk_rows = mpi_size/chunk_cols;

    int chunk_x = mpi_rank % chunk_cols;
    int chunk_y = (mpi_rank / chunk_cols) % chunk_rows;

    int chunk_w = COLS / chunk_cols;
    int chunk_h = ROWS / chunk_rows;

    int neighbor_N = mpi_rank - chunk_cols;
    int neighbor_S = mpi_rank + chunk_cols;
    int neighbor_E = mpi_rank + 1;
    int neighbor_W = mpi_rank - 1;
    
    int local_I_FIX = I_FIX - (chunk_y * chunk_h);
    int local_J_FIX = J_FIX - (chunk_x * chunk_w);

    int contains_fix = (local_I_FIX >= 0 && local_I_FIX < chunk_h)
                       && (local_J_FIX >= 0 && local_J_FIX < chunk_w);

    printf("%d %d %d rank %d size %d  rows %d cols %d  x %d y %d w %d h %d N %d S %d E %d W %d\n",contains_fix,local_I_FIX, local_J_FIX, mpi_rank, mpi_size, chunk_rows, chunk_cols, chunk_x, chunk_y, chunk_w, chunk_h, neighbor_N, neighbor_S, neighbor_E, neighbor_W);

    //MPI_Finalize();
    //return 0;

    int array_w = chunk_w + 2;
    int array_h = chunk_h + 2;

    double* sendbuf_row = (double *) malloc(chunk_w*sizeof(double));
    double* recvbuf_row = (double *) malloc(chunk_w*sizeof(double));
    double* sendbuf_col = (double *) malloc(chunk_h*sizeof(double));
    double* recvbuf_col = (double *) malloc(chunk_h*sizeof(double));

    int tag = 0;

    double **a_old = alloc_matrix(array_w+2, array_h+2); //allocate memory for the matrices + 2 for shared edges
    double **a_new = alloc_matrix(array_w+2, array_h+2);

    init_matrix(a_old, array_w, array_h); //initialize the matrices
    init_matrix(a_new, array_w, array_h);

    if(contains_fix){
      a_old[local_I_FIX+1][local_J_FIX+1] = TEMP;
      a_new[local_I_FIX+1][local_J_FIX+1] = TEMP;
    }

    while (1) {

        //Pass South 
        if(chunk_rows == 1); // no neighbors N or S
          //do nothing
        else if (chunk_y > 0 && chunk_y < (chunk_rows-1)) { //neighbors both N and S
          copy_row_to_buf(a_old, sendbuf_row, chunk_w, chunk_h);
          MPI_Sendrecv(sendbuf_row, chunk_w, MPI_DOUBLE, neighbor_S, tag,
                       recvbuf_row, chunk_w, MPI_DOUBLE, neighbor_N, tag,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          copy_row_from_buf(a_old, recvbuf_row, chunk_w, 0);
        } 
        else if (chunk_y > 0) { //just N neighbor
          MPI_Recv(recvbuf_row, chunk_w, MPI_DOUBLE, neighbor_N, tag,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          copy_row_from_buf(a_old, recvbuf_row, chunk_w, 0);
        } else { //just S neighbor
          copy_row_to_buf(a_old, sendbuf_row, chunk_w, chunk_h);
          MPI_Send(sendbuf_row, chunk_w, MPI_DOUBLE, neighbor_S, tag,
                   MPI_COMM_WORLD);
        }

        //Pass East 
        if(chunk_cols == 1); // no neighbors E or W
          //do nothing
        else if (chunk_x > 0 && chunk_x < (chunk_cols-1)) {
          copy_col_to_buf(a_old, sendbuf_col, chunk_h, chunk_w);
          MPI_Sendrecv(sendbuf_col, chunk_h, MPI_DOUBLE, neighbor_E, tag,
                       recvbuf_col, chunk_h, MPI_DOUBLE, neighbor_W, tag,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          copy_col_from_buf(a_old, recvbuf_col, chunk_h, 0);
        } 
        else if (chunk_x > 0) { //just W neighbor
          MPI_Recv(recvbuf_col, chunk_h, MPI_DOUBLE, neighbor_W, tag,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          copy_col_from_buf(a_old, recvbuf_col, chunk_h, 0);
        } else { //just E neighbor
          copy_col_to_buf(a_old, sendbuf_col, chunk_h, chunk_w);
          MPI_Send(sendbuf_col, chunk_h, MPI_DOUBLE, neighbor_E, tag,
                   MPI_COMM_WORLD);
        }

        //Pass North
        if(chunk_rows == 1); // no neighbors N or S
          //do nothing
        else if (chunk_y > 0 && chunk_y < (chunk_rows-1)) { //neighbors both N and S
          copy_row_to_buf(a_old, sendbuf_row, chunk_w, 1);
          MPI_Sendrecv(sendbuf_row, chunk_w, MPI_DOUBLE, neighbor_N, tag,
                       recvbuf_row, chunk_w, MPI_DOUBLE, neighbor_S, tag,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          copy_row_from_buf(a_old, recvbuf_row, chunk_w, array_h-1);
        } 
        else if (chunk_y > 0) { //just N neighbor
          copy_row_to_buf(a_old, sendbuf_row, chunk_w, 1);
          MPI_Send(sendbuf_row, chunk_w, MPI_DOUBLE, neighbor_N, tag,
                   MPI_COMM_WORLD);
        } else { //just S neighbor
          MPI_Recv(recvbuf_row, chunk_w, MPI_DOUBLE, neighbor_S, tag,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          copy_row_from_buf(a_old, recvbuf_row, chunk_w, array_h-1);
        }
  
        //Pass West 
        if(chunk_cols == 1); // no neighbors E or W
          //do nothing
        else if (chunk_x > 0 && chunk_x < (chunk_cols-1)) {
          copy_col_to_buf(a_old, sendbuf_col, chunk_h, 1);
          MPI_Sendrecv(sendbuf_col, chunk_h, MPI_DOUBLE, neighbor_W, tag,
                       recvbuf_col, chunk_h, MPI_DOUBLE, neighbor_E, tag,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          copy_col_from_buf(a_old, recvbuf_col, chunk_h, array_w-1);
        } 
        else if (chunk_x > 0) { //just W neighbor
          copy_col_to_buf(a_old, sendbuf_col, chunk_h, 1);
          MPI_Send(sendbuf_col, chunk_h, MPI_DOUBLE, neighbor_W, tag,
                   MPI_COMM_WORLD);
        } else { //just E neighbor
          MPI_Recv(recvbuf_col, chunk_h, MPI_DOUBLE, neighbor_E, tag,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          copy_col_from_buf(a_old, recvbuf_col, chunk_h, array_w-1);
        }

        //compute new values and put them into a_new
        compute_new_values(a_old, a_new, array_w, array_h);
        if(contains_fix)
          a_new[local_I_FIX+1][local_J_FIX+1] = TEMP;

        if (DEBUG) {
            printf("a_old is:\n"); //output matrix to screen
            print_matrix(a_old, array_w, array_h);
            printf("a_new is:\n");
            print_matrix(a_new, array_w, array_h);
        }

        //calculate the maximum absolute differences among pairwise
        // differences of old and new matrix elements
        double max_diff = max_abs(a_old, a_new, array_w, array_h);
        int done = (max_diff < EPS);
        int count_done;
        MPI_Allreduce(&done, &count_done, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if(count_done == mpi_size)
            break;

        copy_matrix(a_old, a_new, array_w, array_h); //assign values of a_new to a_old
    }

    double** result_matrix = alloc_matrix(COLS, ROWS);

    int* recvcounts = (int*)malloc(mpi_size*sizeof(int));
    int* displs = (int*)malloc(mpi_size*sizeof(int));
     

    for(int i = 0; i < ROWS; i++)
    {

      if(mpi_rank == root) {
        for(int rank = 0; rank < mpi_size; rank++){
          int rank_x = rank % chunk_cols;
          int rank_y = (rank / chunk_cols) % chunk_rows;
          int rank_i = i - (rank_y * chunk_h);

          if(rank_i >= 0 && rank_i < chunk_h){
            recvcounts[rank] = chunk_w;
            displs[rank] = rank_x*chunk_w;
          } else {
            recvcounts[rank] = 0;
            displs[rank] = 0;
          }
        }
      }
       
      int local_i = i - (chunk_y * chunk_h);
      int sendcount = 0;
      if( local_i >= 0 && local_i < chunk_h){
          copy_row_to_buf(a_new, sendbuf_row, chunk_w, local_i+1);
          sendcount = chunk_w;
      }

      MPI_Gatherv(sendbuf_row, sendcount, MPI_DOUBLE,
                  result_matrix[i], recvcounts, displs, MPI_DOUBLE,
                  root, MPI_COMM_WORLD);
    }

    if (mpi_rank == root) {
      print_matrix_csv(result_matrix, COLS, ROWS);
      /*
      FILE * outfile;
      outfile = fopen ("file.txt", "w+");
       
      printf("\nThe final heat distribution matrix is:\n");
      print_matrix(result_matrix, COLS, ROWS);
      print_matrix_csv(outfile, result_matrix, COLS, ROWS);

      fclose(outfile);
      */
      time(&end);
      double timediff = difftime(end, start);
      printf("execution time (wall): %f s\n", timediff);
    }

    free(sendbuf_row);
    free(recvbuf_row);
    free(sendbuf_col);
    free(recvbuf_col);
    free(a_old);
    free(a_new);
    free(result_matrix);

    MPI_Finalize();
    return 0;
}
