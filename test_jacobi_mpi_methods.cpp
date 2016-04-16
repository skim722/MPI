#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

template <typename T>
void print_vec(vector<T> &foo) {
    cout << "{ ";
    for (int i=0; i < foo.size(); ++i) cout << foo[i] << ", ";
    cout << "}\n";
}

template <typename T>
void print_matrix(vector<T> &foo, int m, int n) {
    for (int i=0; i < m; ++i) {
        cout << "[ ";
        for (int j=0; j < n; ++j) {
            cout << foo[i*n + j] << ", ";
        } cout << "]\n";
    }
}


template <typename T>
void mpi_print_vec(vector<T> &lst, MPI_Comm comm) {
    int num_processors, rank;
    MPI_Comm_size(comm, &num_processors);
    MPI_Comm_rank(comm, &rank);

    MPI_Barrier(comm);
    for (int i=0; i < num_processors; ++i) {
        if (i == rank) {
            cout << "PROCESSOR [" << i << "]:  ";
            print_vec(lst);
            cout << "\n";
        }
        MPI_Barrier(comm);
    }
}

template <typename T>
void mpi_print_vec_cart(vector<T> &lst, MPI_Comm comm) {
    int maxdims=2;
    vector<int> dims(maxdims,0), periods(maxdims,0), coords(maxdims,0);
    MPI_Cart_get(comm, maxdims, &dims[0], &periods[0], &coords[0]);

    if (coords[0] == 0 && coords[1] == 0) {
        cout << "\n\n";
    }

    MPI_Barrier(comm);
    for (int i=0; i < dims[0]; ++i) {
        for (int j=0; j < dims[1]; ++j) {
            if (coords[0] == i && coords[1] == j) {
                cout << "PROCESSOR [" << i << ", " << j << "]:  ";
                print_vec(lst);
                cout << "\n";
            }
            MPI_Barrier(comm);
        }
    }

    if (coords[0] == 0 && coords[1] == 0) {
        cout << "\n\n";
    }
    MPI_Barrier(comm);
}


template <typename T>
void mpi_print_matrix_cart(vector<T> &mat, MPI_Comm comm, int m, int n) {
    int maxdims=2;
    vector<int> dims(maxdims,0), periods(maxdims,0), coords(maxdims,0);
    MPI_Cart_get(comm, maxdims, &dims[0], &periods[0], &coords[0]);

    if (coords[0] == 0 && coords[1] == 0) {
        cout << "\n\n";
    }

    MPI_Barrier(comm);
    for (int i=0; i < dims[0]; ++i) {
        for (int j=0; j < dims[1]; ++j) {
            if (coords[0] == i && coords[1] == j) {
                cout << "PROCESSOR [" << i << ", " << j << "]:\n";
                print_matrix(mat, m, n);
                cout << "\n";
            }
            MPI_Barrier(comm);
        }
    }

    if (coords[0] == 0 && coords[1] == 0) {
        cout << "\n\n";
    }
    MPI_Barrier(comm);
}

vector<int> get_cartesian_coords(MPI_Comm comm, int maxdims=2) {
    vector<int> dims(maxdims,0), periods(maxdims,0), coords(maxdims,0);
    MPI_Cart_get(comm, maxdims, &dims[0], &periods[0], &coords[0]);
    return coords;
}

int get_rowcol_head_rank(MPI_Comm rowcol_comm) {
    int rank0;
    int coords[2] = {0, 0};
    MPI_Cart_rank(rowcol_comm, coords, &rank0);
    return rank0;
}

int get_chunk_size(int n, int q, int i) {
    if (i < n % q) {
        return ceil((double)n/q);
    } else {
        return floor((double)n/q);
    }
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

vector<double> prefilled_matrix(int m, int n) {
    vector<double> A(m * n, 0);
    for (int i=0; i < m; i++) {
        for (int j=0; j < n; ++j) {
            A[i*n + j] = 1.0 + i / 10.0 + j / 100.0;
        }
    }
    return A;
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

void test_mpi_bcast(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int p; MPI_Comm_size(MPI_COMM_WORLD, &p);

    double arr[] = { 213.33, 2355.56, 85.454, 6783.234, -234.67, 33.97, -98789.22, 8798.1, -987098.0021, 232 };
    vector<double> input_vector(arr, arr + sizeof(arr) / sizeof(arr[0]) );
    vector<double> local_vector(4, 0);

    MPI_Scatter(&input_vector[0], local_vector.size(), MPI_DOUBLE, &local_vector[0], local_vector.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    mpi_print_vec(local_vector, MPI_COMM_WORLD);

    MPI_Finalize();
}

void test_distribute_vector(int argc, char *argv[]) {
    // mpic++ test_jacobi_mpi_methods.cpp && mpirun -n 12 a.out

    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int p; MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (p != 12) {
        cerr << "Need 12 processors for this test; exiting!\n";
        MPI_Finalize();
        return;
    }

    double arr[] = { 213.33, 2355.56, 85.454, 6783.234, -234.67, 33.97, -98789.22, 8798.1, -987098.0021, 232 };
    vector<double> input_vector(arr, arr + sizeof(arr) / sizeof(arr[0]) );
    if (rank != 0) {
        input_vector = vector<double>(input_vector.size(), 0);
    }

    if (rank == 2) {
        double arr2[] = { 13.33, 111, 222, 66.234, -234.67, 0.22, -9.22, 8.1, -98.0021, 65 };
        input_vector = vector<double>(arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Create 2D grid comm
    MPI_Comm grid_comm;
    int dims[2] = {3, 4};
    int periods[2] = {0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid_comm);

    // show state of grid
    mpi_print_vec_cart(input_vector, grid_comm);

    // Create column comm
    MPI_Comm column_comm;
    int remain_dims[] = { true, false };
    MPI_Cart_sub(grid_comm, remain_dims, &column_comm);

    // Get chunk sizes per processor and displacements
    vector<int> chunks(dims[0], 0), disps(dims[0], 0);
    for (int i=0; i < dims[0]; ++i) {
        chunks[i] = get_chunk_size(input_vector.size(), dims[0], i);
    }
    for (int i=1; i < chunks.size(); ++i) {
        disps[i] = disps[i-1] + chunks[i-1];
    }

    // Allocate local_vector
    vector<int> coords = get_cartesian_coords(column_comm);
    vector<double> local_vector(get_chunk_size(input_vector.size(), dims[0], coords[0]), 0);

    // ScatterV the elements (not Scatter)
    MPI_Scatterv(&input_vector[0], &chunks[0], &disps[0], MPI_DOUBLE, &local_vector[0], local_vector.size(), MPI_DOUBLE, 0, column_comm);
    mpi_print_vec_cart(local_vector, grid_comm);

    MPI_Comm_free(&column_comm);
    MPI_Comm_free(&grid_comm);
    MPI_Finalize();
}

void test_transpose_bcast_vector(int argc, char *argv[]) {
    // mpic++ test_jacobi_mpi_methods.cpp && mpirun -n 16 a.out
    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int p; MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (p != 16) {
        cerr << "Need 16 processors for this test; exiting!\n";
        MPI_Finalize();
        return;
    }

    double arr[] = { 213.33, 2355.56, 85.454, 6783.234, -234.67, 33.97, -98789.22, 8798.1, -987098.0021, 232 };
    vector<double> input_vector(arr, arr + sizeof(arr) / sizeof(arr[0]) );
    if (rank != 0) {
        input_vector = vector<double>(input_vector.size(), 0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Create 2D grid comm
    MPI_Comm grid_comm;
    int dims[2] = {4, 4};
    int periods[2] = {0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid_comm);

    // show state of grid
    mpi_print_vec_cart(input_vector, grid_comm);

    // Create column comm
    MPI_Comm column_comm;
    int remain_dims[] = { true, false };
    MPI_Cart_sub(grid_comm, remain_dims, &column_comm);

    // Get chunk sizes per processor and displacements
    vector<int> chunks(dims[0], 0), disps(dims[0], 0);
    for (int i=0; i < dims[0]; ++i) {
        chunks[i] = get_chunk_size(input_vector.size(), dims[0], i);
    }
    for (int i=1; i < chunks.size(); ++i) {
        disps[i] = disps[i-1] + chunks[i-1];
    }

    // Allocate local_vector
    vector<int> coords = get_cartesian_coords(column_comm);
    vector<double> local_vector(get_chunk_size(input_vector.size(), dims[0], coords[0]), 0);

    // ScatterV the elements (not Scatter)
    MPI_Scatterv(&input_vector[0], &chunks[0], &disps[0], MPI_DOUBLE, &local_vector[0], local_vector.size(), MPI_DOUBLE, 0, column_comm);
    mpi_print_vec_cart(local_vector, grid_comm);



    // The important part - bcast_vector

    // Get the diagonal (re-calculate the coords in terms of the grid_comm)
    coords = get_cartesian_coords(grid_comm);
    int rank_source, rank_dest;
    MPI_Cart_shift(grid_comm, 1, coords[0], &rank_source, &rank_dest);

    // Copy to diagonal using send() and receive()
    if (coords[1] == 0) {
        MPI_Send(&local_vector[0], local_vector.size(), MPI_DOUBLE, rank_dest, 0, grid_comm);
    } else if (coords[0] == coords[1]) {
        MPI_Status stat;
        MPI_Recv(&local_vector[0], local_vector.size(), MPI_DOUBLE, rank_source, 0, grid_comm, &stat);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_print_vec_cart(local_vector, grid_comm);


    // Copy to vector2, but resize it so that all processors have same vector size to receive
    vector<double> local_vector2(local_vector);
    local_vector2.resize(get_chunk_size(input_vector.size(), dims[0], coords[1]));
    mpi_print_vec_cart(local_vector2, grid_comm);

    // Broadcast along the column, with the diagonal element as the root for each communicator
    MPI_Bcast(&local_vector2[0], local_vector2.size(), MPI_DOUBLE, coords[1], column_comm);
    mpi_print_vec_cart(local_vector2, grid_comm);

    MPI_Comm_free(&column_comm);
    MPI_Comm_free(&grid_comm);
    MPI_Finalize();
}

void test_distribute_matrix_deprecated(int argc, char *argv[]) {
    // mpic++ test_jacobi_mpi_methods.cpp && mpirun -n 4 a.out
    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int p; MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (p != 4) {
        cerr << "Need 4 processors for this test; exiting!\n";
        MPI_Finalize();
        return;
    }

    // Create 2D grid comm and get local processor's cartesian coords
    MPI_Comm grid_comm;
    int dims[2] = {2, 2};
    int periods[2] = {0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid_comm);
    vector<int> coords = get_cartesian_coords(grid_comm);

    // Initialize the full matrix at the root processor
    // int m = 12, n = 7;
    int m = 8, n = 8;
    vector<double> A(m * n, 0);
    if (rank == 0) {
        A = prefilled_matrix(m, n);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_print_matrix_cart(A, grid_comm, m, n);

    // Create the vector MPI datatype for the submatrix of submatrix_m rows * submatrix_n columns
    MPI_Datatype matrix_block_type2, matrix_block_type;
    int submatrix_m = get_chunk_size(m, dims[0], coords[0]);
    int submatrix_n = get_chunk_size(n, dims[1], coords[1]);
    /*
        The vector consists of submatrix_m blocks of submatrix_n consecutive elements,
        with a distance of n elements between the start of each block
    */
    MPI_Type_vector( submatrix_m, submatrix_n, n, MPI_DOUBLE, &matrix_block_type2 );
    MPI_Type_create_resized( matrix_block_type2, 0, sizeof(double), &matrix_block_type);
    MPI_Type_commit(&matrix_block_type);

    int NPROWS = dims[0];
    int NPCOLS = dims[1];

    int disps[NPROWS*NPCOLS];
    for (int i=0; i < NPROWS; i++) {
        for (int j=0; j < NPCOLS; j++) {
            disps[i*NPCOLS + j] = i*n*submatrix_m + j*submatrix_n;
        }
    }

    // Compute the displacements
    vector<int> displacements(dims[0] * dims[1], 0);
    int last_coords[] = { 0, 0 };
    for (int i=0; i < dims[0]; i++) {
        for (int j=0; j < dims[1]; j++) {
            if ((i == 0 && j == 0)) continue;
            int pos = i*dims[1] + j;
            displacements[pos] = displacements[pos - 1] + get_chunk_size(m, dims[0], last_coords[0]);
            last_coords[0] = i;
            last_coords[1] = j;
        }
    }

    mpi_print_vec_cart(displacements, grid_comm);

    /*
    for (int i=0; i < NPROWS; i++) {
        for (int j=0; j < NPCOLS; j++) {
            disps[i*NPCOLS + j] = i*n*submatrix_m + j*submatrix_n;
        }
    }
    */

    vector<int> counts(dims[0] * dims[1], 1);
    vector<double> ALocal(submatrix_m * submatrix_n, 0);
    MPI_Scatterv(&A[0], &counts[0], &displacements[0], matrix_block_type, &ALocal[0], submatrix_m*submatrix_n, MPI_DOUBLE, 0, grid_comm);
    mpi_print_matrix_cart(ALocal, grid_comm, submatrix_m, submatrix_n);

    MPI_Comm_free(&grid_comm);
    MPI_Finalize();
}

void test_distribute_matrix(int argc, char *argv[]) {
    // mpic++ test_jacobi_mpi_methods.cpp && mpirun -n 15 a.out
    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int p; MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (p != 15) {
        cerr << "Need 15 processors for this test; exiting!\n";
        MPI_Finalize();
        return;
    }

    // Create 2D grid comm and get local processor's cartesian coords
    MPI_Comm grid_comm;
    int dims[2] = {3, 5};
    int periods[2] = {0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid_comm);
    vector<int> coords = get_cartesian_coords(grid_comm);

    // Initialize the full matrix at the root processor
    int m = 11, n = 7;
    //int m = 8, n = 8;
    vector<double> A_global(m * n, 0);
    if (rank == 0) {
        A_global = prefilled_matrix(m, n);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_print_matrix_cart(A_global, grid_comm, m, n);


    // Compute the size of the submatrix block to be held by the CURRENT processor
    int submatrix_m = get_chunk_size(m, dims[0], coords[0]);
    int submatrix_n = get_chunk_size(n, dims[1], coords[1]);

    vector<int> submatrix_dims(2, 0); submatrix_dims[0] = submatrix_m; submatrix_dims[1] = submatrix_n;
    mpi_print_vec_cart(submatrix_dims, grid_comm);



    /* PART 1 = Separate the matrix into row-blocks and ScatterV them to the processors in the first column */

    // Create column comm and get coordinates
    MPI_Comm column_comm;
    int remain_dims[] = { true, false };
    MPI_Cart_sub(grid_comm, remain_dims, &column_comm);


    // Compute the sendcounts and displacements
    vector<int> sendcounts(dims[0], 0), displacements(dims[0], 0);
    for (int i=0; i < dims[0]; ++i) {
        sendcounts[i] = get_chunk_size(m, dims[0], i) * n;
        if (i > 0) {
            displacements[i] = displacements[i-1] + sendcounts[i-1];
        }
    }

    vector<double> A_rowblock(submatrix_m * n, 0);
    if (coords[1] == 0) {
        MPI_Scatterv(&A_global[0], &sendcounts[0], &displacements[0], MPI_DOUBLE, &A_rowblock[0], A_rowblock.size(), MPI_DOUBLE, 0, column_comm);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_print_matrix_cart(A_rowblock, grid_comm, submatrix_m, n);



    /* PART 2 = Separate the row-blocks into column subblocks and ScatterV them to the processors in the first row */

    // Create row comm
    MPI_Comm row_comm;
    int remain_dims2[] = { false, true };
    MPI_Cart_sub(grid_comm, remain_dims2, &row_comm);

    /*
        Create an MPI_vector type that represents a COLUMN of the rowblock
        The vector consists of submatrix_m blocks of 1 consecutive elements,
        with a distance of n elements between the start of each block
    */
    MPI_Datatype column_type2, column_type;
    MPI_Type_vector(submatrix_m, 1, n, MPI_DOUBLE, &column_type2);
    MPI_Type_create_resized(column_type2, 0, 1*sizeof(double), &column_type);
    MPI_Type_commit(&column_type);

    // Compute the sendcounts and displacements
    sendcounts.clear(); sendcounts.resize(dims[1], 0);
    displacements.clear(); displacements.resize(dims[1], 0);
    for (int i=0; i < dims[1]; ++i) {
        // Contrary to the first ScatterV, we are sending COLUMNs, not ELEMENTS, so we only need get_chunk_size()
        sendcounts[i] = get_chunk_size(n, dims[1], i);
        if (i > 0) {
            displacements[i] = displacements[i-1] + sendcounts[i-1];
        }
    }

    // Scatterv the column sub-blocks off the rowblock
    vector<double> A_local(submatrix_m * submatrix_n, 0);
    MPI_Scatterv(&A_rowblock[0], &sendcounts[0], &displacements[0], column_type, &A_local[0], A_local.size(), MPI_DOUBLE, 0, row_comm);
    mpi_print_matrix_cart(A_local, grid_comm, submatrix_n, submatrix_m); // Printed as transposed form N x M instead of M x N

    // Transpose it to normal form
    transpose_matrix(A_local, submatrix_n, submatrix_m);
    mpi_print_matrix_cart(A_local, grid_comm, submatrix_m, submatrix_n);

    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&column_comm);
    MPI_Comm_free(&grid_comm);
    MPI_Finalize();
}

int main(int argc, char *argv[]) {
    // test_distribute_vector(argc, argv);
    // test_transpose_bcast_vector(argc, argv);
    // test_distribute_matrix_deprecated(argc, argv);
    test_distribute_matrix(argc, argv);
    return 0;
}

