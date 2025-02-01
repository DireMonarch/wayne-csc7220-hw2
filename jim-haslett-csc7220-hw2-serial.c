/**
 * Copyright 2025 Jim Haslett
 *
 * This work part of a university assignment.  If you are taking the course
 * this work was assigned for, do the right thing, and solve the assignment
 * yourself!
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define OUTFILENAME "result-serial.txt"
#define MATRIX_SIZE 400



double wtime() {
    /**
     * wtime function used to generate a wall time, in seconds, similar to MPI_wtime()
     * 
     * Based entirely on the MPI_wtime implementation:
     * https://github.com/open-mpi/ompi/blob/main/ompi/mpi/c/wtime.c
     */

    double wtime;
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
    wtime  = (double)tp.tv_nsec/1.0e+9;
    wtime += tp.tv_sec;
    return wtime;
}


int main(int argc, char* argv[]) {
    /* matrix and vector */
    int* vector = (int *)malloc(MATRIX_SIZE * sizeof(int));
    int* matrix = (int *)malloc(MATRIX_SIZE * MATRIX_SIZE * sizeof(int));

    /* product vector */
    int* product_vector = (int *)malloc(MATRIX_SIZE * sizeof(int));

    /* Timer variables */
    double start_time, end_time;

    int value;
    for (int i = 0; i < MATRIX_SIZE; i++) {
        value = 1 + i;
        for (int j = 0; j < MATRIX_SIZE; j++) {
            matrix[j*MATRIX_SIZE + i] = value;
        }
        /* initialize vector to all 1 */
        vector[i] = 1;

        /* initialize product_vector to all zeros */
        product_vector[i] = 0;
    }

    /* Record start time */
    start_time = wtime();

    /* Main matrix vector mulitplication */
    for (int j = 0; j < MATRIX_SIZE; j++) {
        for (int i = 0; i < MATRIX_SIZE; i++) {
            product_vector[j] += matrix[j*MATRIX_SIZE + i] * vector[i];
        }
    }

    /* record end time */
    end_time = wtime();

    /* Output result vector to file */
    FILE *outfile; /* File Pointer for output */
    outfile = fopen(OUTFILENAME, "w");
    for (int i = 0; i < MATRIX_SIZE; i++) {
        fprintf(outfile, "%d  ", product_vector[i]);
    }
    fprintf(outfile, "\n");
    fclose(outfile);


    /* calculate times */
    double total_time = end_time - start_time;

    /* Print execution time to screen*/
    printf("\nSerial Matrix-Vector multiply of size %d - time: %f \n\n", MATRIX_SIZE, total_time);

    free(matrix);
    free(vector);
    free(product_vector);

    exit(0);
}
