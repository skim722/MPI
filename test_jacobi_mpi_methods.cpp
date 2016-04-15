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
void mpi_print_vec_cart(vector<T> &lst, MPI_Comm comm, int maxdims=2) {
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


vector<int> get_cartesian_coords(MPI_Comm comm, int maxdims=2) {
    vector<int> dims(maxdims,0), periods(maxdims,0), coords(maxdims,0);
    MPI_Cart_get(comm, maxdims, &dims[0], &periods[0], &coords[0]);
    return coords;
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

void test_distribute_vector(int argc, char *argv[]) {
    // mpic++ test_jacobi_mpi_methods.cpp && mpirun -n 12 a.out

    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int p; MPI_Comm_size(MPI_COMM_WORLD, &p);

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

    // Create 2D 3x4 grid comm
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

void transpose_bcast_vector(int argc, char *argv[]) {
    // mpic++ test_jacobi_mpi_methods.cpp && mpirun -n 16 a.out
    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int p; MPI_Comm_size(MPI_COMM_WORLD, &p);

    double arr[] = { 213.33, 2355.56, 85.454, 6783.234, -234.67, 33.97, -98789.22, 8798.1, -987098.0021, 232 };
    vector<double> input_vector(arr, arr + sizeof(arr) / sizeof(arr[0]) );
    if (rank != 0) {
        input_vector = vector<double>(input_vector.size(), 0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Create 2D 3x4 grid comm
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
    local_vector2.resize(ceil((double)input_vector.size()/dims[0]));
    mpi_print_vec_cart(local_vector2, grid_comm);

    // Broadcast along the column, with the diagonal element as the root for each communicator
    MPI_Bcast(&local_vector2[0], local_vector2.size(), MPI_DOUBLE, coords[1], column_comm);
    mpi_print_vec_cart(local_vector2, grid_comm);

    MPI_Comm_free(&column_comm);
    MPI_Comm_free(&grid_comm);
    MPI_Finalize();
}



int main(int argc, char *argv[]) {
    // test_distribute_vector(argc, argv);
    transpose_bcast_vector(argc, argv);
    return 0;
}

