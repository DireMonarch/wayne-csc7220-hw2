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
#define REPORTFILENAME "report.txt"
#define EXPECTED_PROCESSES 16
#define MESH_DIMS 4
#define LOCAL_MATRIX_SIZE 100


int main(int argc,char* argv[]) {
    /* MPI rank */
    int rank;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);

    /* Fetch rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Fetch number of processes */
    int number_of_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    if(number_of_processes != EXPECTED_PROCESSES){
        printf("Rank %d process expected to see %d process, but saw %d instead.  Exiting!\n", rank, EXPECTED_PROCESSES, number_of_processes);
        MPI_Finalize();
        exit(1);
    }

    /* Create the mesh */
    MPI_Comm comm_mesh;
    int dims[2], periods[2];
    dims[0] = dims[1] = MESH_DIMS;
    periods[0] = periods[1] = 1;
    MPI_Cart_create( MPI_COMM_WORLD , 2 , dims , periods , 1 , &comm_mesh);

    int rank_2d, my_coords[2];
    MPI_Comm_rank( comm_mesh , &rank_2d);
    MPI_Cart_coords( comm_mesh , rank_2d , 2 , my_coords);

    /* Get row and column communicators */
    MPI_Comm comm_row, comm_col;
    int retain_dims[2];
    retain_dims[0] = 0;// This will keep the rows
    retain_dims[1] = 1;
    MPI_Cart_sub(comm_mesh, retain_dims, &comm_row);
    retain_dims[0] = 1;
    retain_dims[1] = 0;
    MPI_Cart_sub(comm_mesh, retain_dims, &comm_col);

    /* Initialize the local matrix and vector */
    int* local_vector = (int *)malloc(LOCAL_MATRIX_SIZE * sizeof(int));
    int* local_matrix = (int *)malloc(LOCAL_MATRIX_SIZE * LOCAL_MATRIX_SIZE * sizeof(int));

    /* Initialize the local product vector */
    int* product_vector = (int *)malloc(LOCAL_MATRIX_SIZE * sizeof(int));

    /* Variables used for communication */
    int other_rank;
    int other_coords[2];

    int value;
    for(int i=0; i<LOCAL_MATRIX_SIZE; i++){
        value = 1 + i + my_coords[1] * LOCAL_MATRIX_SIZE;
        for(int j=0; j<LOCAL_MATRIX_SIZE; j++){
            local_matrix[j*LOCAL_MATRIX_SIZE + i] = value;
        }
        if (my_coords[1] == MESH_DIMS - 1){
            /* only set the vector values to 1 if this is the end column i.e. my_coords[1]==3 */
            local_vector[i] = 1;
        } else {
            /* for all other columns, initialize the vector to 0 */
            local_vector[i] = 0;
        }

        /* set product_vector to all zeros */
        product_vector[i] = 0;
    }

    /* MPI_Recv status variable */
    MPI_Status status;


    /* Communicate the vector parts from the last column to the diagonals */
    /* if current process is in the last column, but not the last row, send 
        vector chunk to the diagonal in the same row */
    if(my_coords[1] == MESH_DIMS-1 && my_coords[0] != my_coords[1]){
        other_coords[0] =  other_coords[1] = my_coords[0];
        MPI_Cart_rank( comm_mesh, other_coords , &other_rank);
        // printf("SEND FROM %d (%d,%d) to %d (%d,%d)\n", rank, my_coords[0], my_coords[1], other_rank, other_coords[0], other_coords[1]);
        MPI_Send( local_vector, LOCAL_MATRIX_SIZE, MPI_INT, other_rank, my_coords[0], comm_mesh);
    }

    if(my_coords[1] != MESH_DIMS-1 && my_coords[0] == my_coords[1]){
        other_coords[0] = my_coords[1];
        other_coords[1] = MESH_DIMS-1;
        MPI_Cart_rank( comm_mesh, other_coords , &other_rank);
        // printf("RECV TO %d (%d,%d) from %d (%d,%d)\n", rank, my_coords[0], my_coords[1], other_rank, other_coords[0], other_coords[1]);
        MPI_Recv( local_vector, LOCAL_MATRIX_SIZE, MPI_INT, other_rank, my_coords[0], comm_mesh, &status);
        // printf("STATUS %d (%d, %d) = %d  %d  %d\n", rank, my_coords[0], my_coords[1], status.MPI_SOURCE, status.MPI_TAG, status.MPI_ERROR);
    } 
    

    /* Communicate the vector parts from the diagonal to the rest of each column */
    other_coords[0] = other_coords[0] = my_coords[1];
    MPI_Cart_rank(comm_col, other_coords, &other_rank);
    MPI_Bcast( local_vector, LOCAL_MATRIX_SIZE, MPI_INT, other_rank, comm_col);


    /* Main matrix vector mulitplication */
    for(int j = 0; j < LOCAL_MATRIX_SIZE; j++){
        for(int i = 0; i < LOCAL_MATRIX_SIZE; i++){
            product_vector[j] += local_matrix[j*LOCAL_MATRIX_SIZE + i] * local_vector[i];
        }
    }


    /* Gather product vector results in last column */
    other_coords[0] = MESH_DIMS - 1;
    MPI_Cart_rank(comm_row, other_coords, &other_rank);
    // int trank;
    // MPI_Comm_rank( comm_row , &trank);
    // printf("My rank %d,  dest %d\n", trank, other_rank);
    int* result = (int *)malloc(LOCAL_MATRIX_SIZE * sizeof(int));
    MPI_Reduce( product_vector, result, LOCAL_MATRIX_SIZE, MPI_INT, MPI_SUM, other_rank, comm_row);

    // printf("rank: %d, %d, %d\n", rank, result[0], product_vector[0]);
    // printf("My 2D Rank: %d  my 2D Coords: (%d, %d), product_vector[0] %d\n", rank_2d, my_coords[0], my_coords[1], result[0]);


    /* collect results */

    /* if we are p0, collect results */
    if(rank == 0){
        free(result);
        result = (int *)malloc(LOCAL_MATRIX_SIZE * MESH_DIMS * sizeof(int));
        other_coords[1] = MESH_DIMS - 1;
        for(int i = 0; i < MESH_DIMS; i++){
            other_coords[0] = i;
            MPI_Cart_rank( comm_mesh, other_coords, &other_rank);
            int chunk = i*LOCAL_MATRIX_SIZE;
            MPI_Recv( result+chunk, LOCAL_MATRIX_SIZE, MPI_INT, other_rank, other_rank, comm_mesh, &status);
            // printf("\n\nChunk: %d  from: %d\n", chunk, other_rank);
            // for(int i = 0; i < LOCAL_MATRIX_SIZE; i++){
            //     printf("%d ", result[i+chunk]);
            // }
        }

        for(int i = 0; i < LOCAL_MATRIX_SIZE*MESH_DIMS; i++){
            printf("%d ", result[i]);
        }
    }

    /* if we are in the last column, send results to p0 */
    if(my_coords[1] == MESH_DIMS-1){
        // if(rank == 15){
        //     printf("\n\nSending from %d\n", rank);
        //     for(int i = 0; i < LOCAL_MATRIX_SIZE; i++){
        //         printf("%d ", result[i]);
        //     }
        // }
        MPI_Send( result, LOCAL_MATRIX_SIZE, MPI_INT, 0, rank, comm_mesh);
    }

    free(local_matrix);
    free(local_vector);
    free(product_vector);
    free(result);

    MPI_Finalize();
    exit(0);
}