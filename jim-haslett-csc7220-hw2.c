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
#include "mpi.h"

#define OUTFILENAME "result.txt"
#define EXPECTED_PROCESSES 16
#define MESH_DIMS 4
#define LOCAL_MATRIX_SIZE 100


int main(int argc, char* argv[]) {
    /* MPI rank */
    int rank;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);

    /* Fetch rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Fetch number of processes */
    int number_of_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    /* If number of processes is not the number expected, quit with an error */
    if (number_of_processes != EXPECTED_PROCESSES) {
        printf("Rank %d process expected to see %d process, but saw %d instead.  Exiting!\n", rank, EXPECTED_PROCESSES, number_of_processes);
        MPI_Finalize();
        exit(1);
    }

    /* Create the mesh communicator */
    MPI_Comm comm_mesh;  // variable used to hold the mesh communicator
    int dims[2], periods[2];
    dims[0] = dims[1] = MESH_DIMS;  // number of processes in each dimension
    periods[0] = periods[1] = 1;  // make the mesh wraparound
    MPI_Cart_create(MPI_COMM_WORLD , 2 , dims , periods , 1 , &comm_mesh);

    /* fetch this processes rank and coordinates in the new mesh communicator */
    int rank_2d, my_coords[2];
    MPI_Comm_rank(comm_mesh , &rank_2d);
    MPI_Cart_coords(comm_mesh , rank_2d , 2 , my_coords);

    /* Get row and column communicators */
    MPI_Comm comm_row, comm_col;
    int retain_dims[2];  // used to tell MPI which dimensions to keep
    /**
     * This will create the row communcators.
     * It is somewhat counterintuitive, in that we tell MPI_Cart_sub which dimensions
     * to *keep* not which dimension we want to build communicators from.  In this case
     * dimension 0, or the rows, we set to 0 (false) and dimension 1, or columns,
     * we set to 1 (true).  MPI will *keep* the column dimensions, and separate
     * on the rows, creating new communicators for each row.
     */
    retain_dims[0] = 0;
    retain_dims[1] = 1;
    MPI_Cart_sub(comm_mesh, retain_dims, &comm_row);

    /* This will create the row communcators. */
    retain_dims[0] = 1;
    retain_dims[1] = 0;
    MPI_Cart_sub(comm_mesh, retain_dims, &comm_col);

    /* Create the local matrix and vector - Inputs */
    int* local_vector = (int *)malloc(LOCAL_MATRIX_SIZE * sizeof(int));
    int* local_matrix = (int *)malloc(LOCAL_MATRIX_SIZE * LOCAL_MATRIX_SIZE * sizeof(int));

    /* Create the local product vector - Output */
    int* product_vector = (int *)malloc(LOCAL_MATRIX_SIZE * sizeof(int));

    /* Variables used for MPI communication */
    int other_rank;
    int other_coords[2];

    /* Timer variables */
    double start_time, end_time;

    int value;  // temporary variable to hold each column calculation.
    for (int i = 0; i < LOCAL_MATRIX_SIZE; i++) {
        /**
         * Per the assignment, the ith column in each row of the matrix has the
         * same value, therefore, we only need to calculate that value once, and
         * set it on each row.  This is a column first loop.
         */
        value = 1 + i + my_coords[1] * LOCAL_MATRIX_SIZE;
        for (int j = 0; j < LOCAL_MATRIX_SIZE; j++) {
            local_matrix[j*LOCAL_MATRIX_SIZE + i] = value;
        }
        /**
         * Since only the last column (MESH_DIMS - 1) contains the input vector,
         * we are setting the local_vector variable to 1, only if this process is
         * in the last column.  Otherwise we initialize the value to 0.
         */
        if (my_coords[1] == MESH_DIMS - 1) {
            /* only set the vector values to 1 if this is the end column i.e. my_coords[1] == 3 */
            local_vector[i] = 1;
        } else {
            /* for all other columns, initialize the vector to 0 */
            local_vector[i] = 0;
        }

        /* initialize product_vector to all zeros */
        product_vector[i] = 0;
    }

    /* MPI_Recv status variable */
    MPI_Status status;

    /* Record start time */
    start_time = MPI_Wtime();

    /**
     * Communicate the vector parts from the last column to the diagonals.
     * If the current process is in the last column, but not the last row, send 
     * vector chunk to the diagonal in the same row.  No need to send the last
     * vector chunk to itself, it already has it!
     */
    if (my_coords[1] == MESH_DIMS-1 && my_coords[0] != my_coords[1]) {
        other_coords[0] =  other_coords[1] = my_coords[0];
        MPI_Cart_rank(comm_mesh, other_coords , &other_rank);
        MPI_Send(local_vector, LOCAL_MATRIX_SIZE, MPI_INT, other_rank, my_coords[0], comm_mesh);
    }

    /* The receive end of the send above. */
    if (my_coords[1] != MESH_DIMS-1 && my_coords[0] == my_coords[1]) {
        other_coords[0] = my_coords[1];
        other_coords[1] = MESH_DIMS-1;
        MPI_Cart_rank(comm_mesh, other_coords , &other_rank);
        MPI_Recv(local_vector, LOCAL_MATRIX_SIZE, MPI_INT, other_rank, my_coords[0], comm_mesh, &status);
    }


    /* Communicate the vector parts from the diagonal to the rest of each column */
    other_coords[0] = other_coords[0] = my_coords[1];
    MPI_Cart_rank(comm_col, other_coords, &other_rank);
    MPI_Bcast(local_vector, LOCAL_MATRIX_SIZE, MPI_INT, other_rank, comm_col);


    /* Main matrix vector mulitplication */
    for (int j = 0; j < LOCAL_MATRIX_SIZE; j++) {
        for (int i = 0; i < LOCAL_MATRIX_SIZE; i++) {
            product_vector[j] += local_matrix[j*LOCAL_MATRIX_SIZE + i] * local_vector[i];
        }
    }

    /* Reduce (sum) product vector results to the last column */
    other_coords[0] = MESH_DIMS - 1;
    MPI_Cart_rank(comm_row, other_coords, &other_rank);
    /* create result array to hold row results */
    int* result = (int *)malloc(LOCAL_MATRIX_SIZE * sizeof(int));
    MPI_Reduce(product_vector, result, LOCAL_MATRIX_SIZE, MPI_INT, MPI_SUM, other_rank, comm_row);

    /* collect results  into processo 0*/

    /* if we are in the last column, send results to p0 */
    if (my_coords[1] == MESH_DIMS-1) {
        MPI_Send(result, LOCAL_MATRIX_SIZE, MPI_INT, 0, rank, comm_mesh);
    }

    /* if we are p0, collect results */
    if (rank == 0) {
        /**
         * result variable (int*) is not used by p0, free any allocated memory
         * and allocate a block to hold the combined result.
         */
        free(result);
        result = (int *)malloc(LOCAL_MATRIX_SIZE * MESH_DIMS * sizeof(int));
        other_coords[1] = MESH_DIMS - 1;  // p0 is alwasy collecting from last column
        /* loop through each row, and collect the results */
        for (int i = 0; i < MESH_DIMS; i++) {
            other_coords[0] = i;
            MPI_Cart_rank(comm_mesh, other_coords, &other_rank);
            int chunk = i*LOCAL_MATRIX_SIZE;  // offset for location to drop results into.
            MPI_Recv(result+chunk, LOCAL_MATRIX_SIZE, MPI_INT, other_rank, other_rank, comm_mesh, &status);
        }

        /* record end time */
        end_time = MPI_Wtime();


        /* calculate times */
        double total_time = end_time - start_time;

        /* Output result vector to file */
        FILE *outfile; /* File Pointer for output */
        outfile = fopen(OUTFILENAME, "w");
        for (int i = 0; i < LOCAL_MATRIX_SIZE*MESH_DIMS; i++) {
            fprintf(outfile, "%d  ", result[i]);
        }
        fprintf(outfile, "\n");
        fclose(outfile);


        /* Print execution time to screen*/
        printf("\nMPI Matrix-Vector multiply of size %d, %d processes - time: %f \n\n", LOCAL_MATRIX_SIZE*MESH_DIMS, number_of_processes, total_time);
    }


    /* Free allocated memory */
    free(local_matrix);
    free(local_vector);
    free(product_vector);
    free(result);

    /* finalize MPI and exit */
    MPI_Finalize();
    exit(0);
}
