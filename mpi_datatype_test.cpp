/*
    Test program to check that MPI_Datatype for MyStruct works.  To run:

        mpic++ test.cpp
        mpirun -n 1 a.out
*/

#include <mpi.h>
#include <iostream>

using namespace std;

typedef struct MyStruct {
    unsigned int key;
    double d;
    char e[4];
} MyStruct;

MPI_Datatype mystruct_get_mpi_type() {
    // use MPI commands to create a custom data type for MyStruct
    MPI_Datatype type, tmp_type;

    // TODO: create the MPI datatype for MyStruct
    MyStruct x;
    MPI_Aint base, adr_key, adr_d, adr_e;
    MPI_Get_address(&x, &base);
    MPI_Get_address(&x.key, &adr_key);
    MPI_Get_address(&x.d, &adr_d);
    MPI_Get_address(&x.e[0], &adr_e);

    MPI_Aint     disps[3] = {adr_key - base, adr_d - base, adr_e - base};
    MPI_Aint     extent   = sizeof(x);
    int          blens[3] = {1, 1, 4};
    MPI_Datatype types[3] = {MPI_UNSIGNED, MPI_DOUBLE, MPI_CHAR};

    MPI_Type_create_struct(3, blens, disps, types, &tmp_type);
    MPI_Type_create_resized(tmp_type, 0, extent, &type);

    MPI_Type_commit(&type);
    return type;
}

void print_MyStruct_arr(MyStruct *arr, size_t arr_size) {
    for (int i=0; i < arr_size; ++i) {
        const MyStruct &x = arr[i];
        cout << "{ " << x.key << ", " << x.d << ", '";
        for (int j=0; j < 4; ++j) {
            cout << x.e[j];
        } cout << "'}\n";
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        MyStruct foos[] = {
            { 42, 3.14, "xyz" },
            { 5787, 4.49, "foo" },
            { 9873, -3.2, "zip" },
        };
        MyStruct bars[3];

        MPI_Datatype type = mystruct_get_mpi_type();

        cout << "Before MPI_Sendrecv call:\n";
        print_MyStruct_arr(bars, 3);

        MPI_Sendrecv(&foos, 3, type, 0, 0, &bars, 3, type, 0, 0, MPI_COMM_SELF, MPI_STATUS_IGNORE);

        cout << "After MPI_Sendrecv call:\n";
        print_MyStruct_arr(bars, 3);
    }

    MPI_Finalize();
    return 0;
}