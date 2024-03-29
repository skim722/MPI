/**
 * @file    mpi_jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements MPI functions for distributing vectors and matrixes,
 *          parallel distributed matrix-vector multiplication and Jacobi's
 *          method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "mpi_jacobi.h"
#include "jacobi.h"
#include "utils.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>

using namespace std;

/*
 * TODO: Implement your solutions here
 */

template <typename T>
void print_matrix(T *foo, int m, int n) {
    cout << "\n";
    for (int i=0; i < m; ++i) {
        cout << "[ ";
        for (int j=0; j < n; ++j) {
            cout << foo[i*n + j] << ", ";
        } cout << "]\n";
    } cout << "\n\n";
}

template <typename T>
void print_vec(T *foo, int n) {
    cout << "\n{ ";
    for (int i=0; i < n; ++i) cout << foo[i] << ", ";
    cout << "}\n\n\n";
}

template <typename T>
void print_vec(vector<T> &foo) {
    print_vec(&foo[0], foo.size());
}

template <typename T>
void print_matrix(vector<T> &foo, int m, int n) {
    print_matrix(&foo[0], m, n);
}

vector<int> get_cartesian_coords(MPI_Comm comm, int maxdims) {
    int rank; MPI_Comm_rank(comm, &rank);
    vector<int> coords(maxdims, 0);
    MPI_Cart_coords(comm, rank, maxdims, &coords[0]);
    return coords;
}

vector<int> get_comm_dims(MPI_Comm comm, int maxdims=2) {
    vector<int> dims(maxdims,0), periods(maxdims,0), coords(maxdims,0);
    MPI_Cart_get(comm, maxdims, &dims[0], &periods[0], &coords[0]);
    return dims;
}

void get_submatrix_size(MPI_Comm comm, int n, int *row, int *col) {
    vector<int> dims(2,0), periods(2,0), coords(2,0);
    MPI_Cart_get(comm, 2, &dims[0], &periods[0], &coords[0]);

    if (coords[0] < n % dims[0]) {
        *row = ceil((double)n/dims[0]);
    } else {
        *row = floor((double)n/dims[0]);
    }

    if (coords[1] < n % dims[1]) {
        *col = ceil((double)n/dims[1]);
    } else {
        *col = floor((double)n/dims[1]);
    }
}

int get_chunk_size(int n, int q, int i) {
    if (i < n % q) {
        return ceil((double)n/q);
    } else {
        return floor((double)n/q);
    }
}

template <typename T>
void transpose_matrix(vector<T> &matrix, int current_m, int current_n) {
    vector<T> matrix2(matrix.size(), 0);
    for (int i = 0; i < current_m; ++i) {
       for (int j = 0; j < current_n; ++j) {
            matrix2[j*current_m + i] = matrix[i*current_n + j];
        }
    }
    matrix = matrix2;
}

void distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm) {
    // Get MPI info - grid dimensions and cartesian coords of current processor
    int maxdims = 2;
    vector<int> dims(maxdims,0), periods(maxdims,0), coords(maxdims,0);
    MPI_Cart_get(comm, maxdims, &dims[0], &periods[0], &coords[0]);

    // Create the column communicator
    MPI_Comm column_comm;
    int remain_dims[] = { true, false };
    MPI_Cart_sub(comm, remain_dims, &column_comm);

    // Determine the send_counts and displacements to each processor in the column
    vector<int> send_counts(dims[0], 0), send_displacements(dims[0], 0);
    for (int i=0; i < dims[0]; ++i) {
        send_counts[i] = get_chunk_size(n, dims[0], i);
    }
    for (int i=1; i < send_counts.size(); ++i) {
        send_displacements[i] = send_displacements[i-1] + send_counts[i-1];
    }

    // Allocate local_vector (use get_chunk_size to determine the size needed for the current processor)
    // Using coords without looking it up a second time using column_comm is fine because
    // in both 2D and 1D communicator, coords[0] is the same
    int local_vector_len = get_chunk_size(n, dims[0], coords[0]);
    *local_vector = new double[ local_vector_len ];

    // ScatterV the elements (not Scatter) ONLY on the first column of processors, since input_vector is a valid pointer only in processor (0,0)
    if (coords[1] == 0) {
        MPI_Scatterv(input_vector, &send_counts[0], &send_displacements[0], MPI_DOUBLE, *local_vector, local_vector_len, MPI_DOUBLE, 0, column_comm);
    }

    MPI_Comm_free(&column_comm);
}

void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm) {
    // Get MPI info - grid dimensions and cartesian coords of current processor
    int maxdims = 2;
    vector<int> dims(maxdims,0), periods(maxdims,0), coords(maxdims,0);
    MPI_Cart_get(comm, maxdims, &dims[0], &periods[0], &coords[0]);

    // Create the column communicator
    MPI_Comm column_comm;
    int remain_dims[] = { true, false };
    MPI_Cart_sub(comm, remain_dims, &column_comm);

    // Determine the send_counts and displacements to each processor in the column
    vector<int> receive_counts(dims[0], 0), receive_displacements(dims[0], 0);
    for (int i=0; i < dims[0]; ++i) {
        receive_counts[i] = get_chunk_size(n, dims[0], i);
    }
    for (int i=1; i < receive_counts.size(); ++i) {
        receive_displacements[i] = receive_displacements[i-1] + receive_counts[i-1];
    }

    // Compute local_vector's length (use get_chunk_size to determine the size needed for the current processor)
    int local_vector_len = get_chunk_size(n, dims[0], coords[0]);

    // GatherV the elements (not Gather) ONLY on the first column of processors, since input_vector is a valid pointer only in processor (0,0)
    if (coords[1] == 0) {
        MPI_Gatherv(local_vector, local_vector_len, MPI_DOUBLE, output_vector, &receive_counts[0], &receive_displacements[0], MPI_DOUBLE, 0, column_comm);
    }

    MPI_Comm_free(&column_comm);
}

void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm) {
    // Get MPI info - grid dimensions and cartesian coords of current processor
    int maxdims = 2;
    vector<int> dims(maxdims,0), periods(maxdims,0), coords(maxdims,0);
    MPI_Cart_get(comm, maxdims, &dims[0], &periods[0], &coords[0]);

    // Since we are provided an N x N matrix, M = N, but we refer to M for consistency
    int m = n;
    int submatrix_m = get_chunk_size(m, dims[0], coords[0]);
    int submatrix_n = get_chunk_size(n, dims[1], coords[1]);

    // Create the column communicator
    MPI_Comm column_comm;
    int remain_dims[] = { true, false };
    MPI_Cart_sub(comm, remain_dims, &column_comm);

    // Create the row communicator
    MPI_Comm row_comm;
    int remain_dims2[] = { false, true };
    MPI_Cart_sub(comm, remain_dims2, &row_comm);


    /*
        PART 1 = Separate the matrix into row-blocks and ScatterV them to the processors in the first column.

        BEFORE: Processor (0,0) contains the global matrix:

            (0,0):
              ________________________________________  _
             |                                        | |  ceil(n/q)
             |                                        | |
             |----------------------------------------| -
             |                                        | ...
             |                                        |
             |----------------------------------------| -
             |                                        | |  floor(n/q)
             |                                        | |
             |----------------------------------------| -

        AFTER: Processor (i,0) contains row-blocks of the global matrix:

            (0,0):
              ________________________________________  _
             |                                        | |  ceil(n/q)
             |                                        | |
             |----------------------------------------| -

            (1,0):
             |----------------------------------------| -
             |                                        | ...
             |                                        |
             |----------------------------------------| -

            (q,0):
             |----------------------------------------| -
             |                                        | |  floor(n/q)
             |                                        | |
             |----------------------------------------| -
    */

    // Compute the sendcounts and displacements for the row-blocks
    vector<int> sendcounts(dims[0], 0), displacements(dims[0], 0);
    for (int i=0; i < dims[0]; ++i) {
        sendcounts[i] = get_chunk_size(m, dims[0], i) * n;
        if (i > 0) {
            displacements[i] = displacements[i-1] + sendcounts[i-1];
        }
    }

    // ScatterV the row-blocks (only the first column of processors participate)
    vector<double> A_rowblock(submatrix_m * n, 0);
    if (coords[1] == 0) {
        MPI_Scatterv(input_matrix, &sendcounts[0], &displacements[0], MPI_DOUBLE, &A_rowblock[0], A_rowblock.size(), MPI_DOUBLE, 0, column_comm);
    }

    // Need to synchronize here since some processors will not be participating in the first ScatterV
    MPI_Barrier(MPI_COMM_WORLD);


    /*
        PART 2 = IN PARALLEL FOR EACH ROW, separate the row-blocks into column subblocks and ScatterV them to the processors in the same row

        BEFORE: Processor (i,0) contains the a row-block:

            (i,0):
             ________________________________________
            |            |             |      |      |
            |            |             |      |      |
            |------------|-------------|------|------|

        AFTER: Processor (i,j) contains the a column-subblock (the processor's block of the full matrix):

            (i,0):          (i,1):           (i,q-1):   (i,q):
             ____________    _____________    ______    ______
            |            |  |             |  |      |  |      |
            |            |  |             |  |      |  |      |
            |------------|  |-------------|  |------|  |------|

        HOWEVER, because we will be using a COLUMN MPI_vector_type to perform the ScatterV here, the elements will appear in transposed form!
        Hence we need to perform matrix transpose in the end!
    */

    /*
        Create an MPI_vector type that represents a COLUMN of the row-block
        The vector consists of submatrix_m blocks of 1 consecutive elements,
        with a distance of n elements between the start of each block
    */
    MPI_Datatype column_type2, column_type;
    MPI_Type_vector(submatrix_m, 1, n, MPI_DOUBLE, &column_type2);
    MPI_Type_create_resized(column_type2, 0, 1*sizeof(double), &column_type);
    MPI_Type_commit(&column_type);

    // Compute the sendcounts and displacements for the column-subblocks
    sendcounts.clear(); sendcounts.resize(dims[1], 0);
    displacements.clear(); displacements.resize(dims[1], 0);
    for (int i=0; i < dims[1]; ++i) {
        // Contrary to the first ScatterV, we are sending COLUMNs, not ELEMENTS, so we only need get_chunk_size() instead of get_chunk_size()*n
        sendcounts[i] = get_chunk_size(n, dims[1], i);
        if (i > 0) {
            displacements[i] = displacements[i-1] + sendcounts[i-1];
        }
    }

    // ScatterV the column-subblocks off the row-block
    vector<double> A_local(submatrix_m * submatrix_n, 0);
    MPI_Scatterv(&A_rowblock[0], &sendcounts[0], &displacements[0], column_type, &A_local[0], A_local.size(), MPI_DOUBLE, 0, row_comm);

    // Allocate local_matrix and fill it with A_local-transpose
    *local_matrix = new double[ A_local.size() ];
    for (int i = 0; i < submatrix_n; ++i) {
       for (int j = 0; j < submatrix_m; ++j) {
            (*local_matrix)[j*submatrix_n + i] = A_local[i*submatrix_m + j];
        }
    }

    // Need to synchronize here since processors of different-sized matrices for transposing
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&column_comm);
}

void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm) {
    // Get MPI info - grid dimensions and cartesian coords of current processor
    int maxdims = 2;
    vector<int> dims(maxdims,0), periods(maxdims,0), coords(maxdims,0);
    MPI_Cart_get(comm, maxdims, &dims[0], &periods[0], &coords[0]);

    // Get the diagonal processor - for a processor with coordinate (k, 0), the diagonal element is k distance away
    int rank_source, rank_dest;
    MPI_Cart_shift(comm, 1, coords[0], &rank_source, &rank_dest);
    int col_vector_len = get_chunk_size(n, dims[0], coords[0]);

    // Send data from the first column to the respective diagonal processors
    if (coords[1] == 0) {
        // If this processor is on the first column (i.e. coordinate (k, 0)), then send()
        MPI_Send(col_vector, col_vector_len, MPI_DOUBLE, rank_dest, 0, comm);
    }

    // Do not use an else-if, because the processor with rank (0,0) needs to send to itself
    if (coords[0] == coords[1]) {
        /*
            If this processor is on the diagonal (i.e. coordinate (k, k)), then receive()
            Copy to the ROW vector since the COL pointer is invalid for processors not in the first column
        */
        MPI_Status stat;
        MPI_Recv(row_vector, col_vector_len, MPI_DOUBLE, rank_source, 0, comm, &stat);
    }

    // Because not all processors participate, we need a barrier here
    MPI_Barrier(MPI_COMM_WORLD);

    // Create the column communicator
    MPI_Comm column_comm;
    int remain_dims[] = { true, false };
    MPI_Cart_sub(comm, remain_dims, &column_comm);

    /*
        Processors in the same column need to know the size of the data being sent
        by the diagonal processor of the column.  Because processor (i,j) is on
        column j, the diagonal processor (the root) is the jth element in the column
        communicator (or row j in the grid communicator), so we use coords[1] instead of coords[0]
    */
    int row_vector_len = get_chunk_size(n, dims[0], coords[1]);

    /*
        Broadcast along the column, with the diagonal element as the root for each communicator
        We use coords[1] to indicate the root processor for the broadcast in the particular column.
        Processor (i,j) is on column j, and so the diagonal processor (the root) will
        be the jth processor in the column communicator
    */
    MPI_Bcast(row_vector, row_vector_len, MPI_DOUBLE, coords[1], column_comm);

    // Deallocate the column communicator
    MPI_Comm_free(&column_comm);
}

void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm) {
    // Get MPI info - grid dimensions and cartesian coords of current processor
    int maxdims = 2;
    vector<int> dims(maxdims,0), periods(maxdims,0), coords(maxdims,0);
    MPI_Cart_get(comm, maxdims, &dims[0], &periods[0], &coords[0]);

    // Create the row communicator
    MPI_Comm row_comm;
    int remain_dims2[] = { false, true };
    MPI_Cart_sub(comm, remain_dims2, &row_comm);

    // Since we are provided an N x N matrix, M = N, but we refer to M for consistency
    int m = n;
    int submatrix_m = get_chunk_size(m, dims[0], coords[0]);
    int submatrix_n = get_chunk_size(n, dims[1], coords[1]);

    /*
        Distribute local_x over the other processors using transpose_bcast_vector()
        Because local_x is only valid on processors in the first column, we need a
        temporary place to store it for all the other processors when calling transpose_bcast_vector()
        The size of each processor's submatrix is submatrix_m x submatrix_n,
        so the length of the process's portion of x (local_x_chunked, NOT local_x)
        should be submatrix_n.  Likewise, the length of the process's portion of y (local_y_chunked, NOT local_y)
        should be submatrix_m.

        NOTE that local_x and local_y refer to the portions of x and y in the processors on the FIRST column of the grid,
        while local_x_chunked and local_y_chunked refer to the portions of x and y in EACH of the processors in the grid

        Based on the implementation of mpi_matrix_vector_mult(),
        it's probably safe to use local_y as well, but we are just being explicit here
    */
    vector<double> local_x_chunked(submatrix_n, 0), local_y_chunked(submatrix_m, 0);
    transpose_bcast_vector(n, local_x, &local_x_chunked[0], comm);

    /*
        Perform the local matrix multiplication and save the local chunk of y to local_y_chunked
    */
    matrix_vector_mult(submatrix_m, submatrix_n, local_A, &local_x_chunked[0], &local_y_chunked[0]);

    /*
        Along each row, the processors will all have a local_y_chunked of length submatrix_m elements
        that needs to be bulk-reduced into the final local_y that is in processor (i,0) (rank 0 in
        the row communicator)
    */
    MPI_Reduce(&local_y_chunked[0], local_y, submatrix_m, MPI_DOUBLE, MPI_SUM, 0, row_comm);
    MPI_Comm_free(&row_comm);
}

void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x, MPI_Comm comm, int max_iter, double l2_termination) {
    // Solves Ax = b using the iterative jacobi method
    // Get MPI info - grid dimensions and cartesian coords of current processor
    int maxdims = 2;
    vector<int> dims(maxdims,0), periods(maxdims,0), coords(maxdims,0);
    MPI_Cart_get(comm, maxdims, &dims[0], &periods[0], &coords[0]);

    // Get the diagonal processor - for a processor with coordinate (k, k), the diagonal element is k distance away to the LEFT
    int rank_source, rank_dest;
    MPI_Cart_shift(comm, 1, -coords[0], &rank_source, &rank_dest);

    // Since we are provided an N x N matrix, M = N, but we refer to M for consistency
    int m = n;
    int submatrix_m = get_chunk_size(m, dims[0], coords[0]);
    int submatrix_n = get_chunk_size(n, dims[1], coords[1]);


    // Initialize D and R - allocate D and copy A to R
    vector<double> local_R(local_A, local_A + (submatrix_m*submatrix_n)), local_D(submatrix_m, 0);
    for (int i=0; i < submatrix_m; ++i) {
        if (coords[0] == coords[1]) {
            // Populate the diagonal
            local_D[i] = local_A[i * submatrix_m + i];
            // Set diagonal of R to be zero
            local_R[i * submatrix_m +i] = 0.0;
        }
        // Initialze x to be zero
        local_x[i] = 0;
    }


    // Distribute D from the diagonal processors to the first column processors
    if (coords[0] == coords[1]) {
        // If this processor is on the diagonal (i.e. coordinate (k, k)), then send()
        MPI_Send(&local_D[0], local_D.size(), MPI_DOUBLE, rank_dest, 0, comm);
    }

    // Do not use an else-if, because the processor with rank (0,0) needs to send to itself
    if (coords[1] == 0) {
        // Else if this processor is on the first column (i.e. coordinate (k, 0)), then receive()
        MPI_Status stat;
        MPI_Recv(&local_D[0], local_D.size(), MPI_DOUBLE, rank_source, 0, comm, &stat);
    }

    // Because not all processors participate, we need a barrier here
    MPI_Barrier(MPI_COMM_WORLD);


    double l2_termination_squared = pow(l2_termination, 2);
    vector<double> local_Rx(submatrix_m, 0), local_Ax(submatrix_m, 0);
    for (int iter=0; iter < max_iter; ++iter) {
        /*
            Compute Rx = R * x
            After this step, only processors of the first column will have a valid local_Rx
        */
        distributed_matrix_vector_mult(n, &local_R[0], local_x, &local_Rx[0], comm);

        /*
            Compute x = (b - rx) / D; this is purely local on the first column
            Processors on other columns will also run this, but this is fine, because
            local_x has been allocated for all processors - the local_x computed in the other columns
            will simply be discarded/ignored
        */
        for (int i=0; i < submatrix_m; ++i) {
            local_x[i] = (local_b[i] - local_Rx[i]) / local_D[i];
        }

        /*
            Compute Ax = A * x
            After this step, only processors of the first column will have a valid local_Ax
        */
        distributed_matrix_vector_mult(n, local_A, local_x, &local_Ax[0], comm);

        /*
            Compute || b - Ax ||  but only if the processor is in the first column
            We need to limit computation to the first column, because we will be
            performing an Allreduce next and we don't want to pollute local_l2_norm_squared
        */
        double local_l2_norm_squared = 0;
        if (coords[1] == 0) {
            for (int i=0; i < submatrix_m; ++i) {
                local_l2_norm_squared += pow(local_b[i] - local_Ax[i], 2);
            }
        }

        /*
            Sum the local_l2_norm_squared values of all processors
            An Allreduce is needed so that the value is broadcast to all processors,
            so that all processors will know whether to stop or continue with the loop
        */
        double global_l2_norm_squared = 0;
        MPI_Allreduce(&local_l2_norm_squared, &global_l2_norm_squared, 1, MPI_DOUBLE, MPI_SUM, comm);

        /*
            If we observe convergence, then we prematurely terminate the algorithm
            We compare to the square of the l2 norm instead of the l2 norm itself to avoid extra computation
        */
        if (global_l2_norm_squared <= l2_termination_squared) return;
    }
}


// wraps the distributed matrix vector multiplication
void mpi_matrix_vector_mult(const int n, double* A,
                            double* x, double* y, MPI_Comm comm)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_x = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);

    distribute_vector(n, &x[0], &local_x, comm);

    // allocate local result space
    double* local_y = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);

    // gather results back to rank 0
    gather_vector(n, local_y, y, comm);
}

// wraps the distributed jacobi function
void mpi_jacobi(const int n, double* A, double* b, double* x, MPI_Comm comm,
                int max_iter, double l2_termination)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_b = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &b[0], &local_b, comm);

    // allocate local result space
    double* local_x = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_jacobi(n, local_A, local_b, local_x, comm, max_iter, l2_termination);

    // gather results back to rank 0
    gather_vector(n, local_x, x, comm);
}
