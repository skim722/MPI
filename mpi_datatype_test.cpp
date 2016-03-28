/*
    Test program to check that MPI_Datatype for MyStruct works.  To run:

        mpi_datatype_test.cpp
        mpirun -n 1 a.out
*/

#include <mpi.h>
#include <iostream>

using namespace std;

struct MyStruct {
    unsigned int key;
    double d;
    char e[4];
};

struct MyStructOpt {
    // Order the struct member variables as such to remove 8 bytes of padding (reduce struct offset (sizeof()) from 24 to 16 bytes)
    double d;           // 8 bytes
    unsigned int key;   // 4 bytes
    char e[4];          // 4 bytes
};


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

MPI_Datatype mystruct_opt_get_mpi_type() {
    // use MPI commands to create a custom data type for MyStructOpt
    MPI_Datatype type, tmp_type;

    // TODO: create the MPI datatype for MyStructOpt
    MyStructOpt x;
    MPI_Aint base, adr_d, adr_key, adr_e;
    MPI_Get_address(&x, &base);
    MPI_Get_address(&x.d, &adr_d);
    MPI_Get_address(&x.key, &adr_key);
    MPI_Get_address(&x.e[0], &adr_e);

    MPI_Aint     disps[3] = {adr_d - base, adr_key - base, adr_e - base};
    MPI_Aint     extent   = sizeof(x);
    int          blens[3] = {1, 1, 4};
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_UNSIGNED, MPI_CHAR};

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

void print_MyStructOpt_arr(MyStructOpt *arr, size_t arr_size) {
    for (int i=0; i < arr_size; ++i) {
        const MyStructOpt &x = arr[i];
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
        {
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

            int size; MPI_Type_size(type, &size);
            cout << "Size of MPI_Datatype for MyStruct:  " << size << "\n";
        }

        {
            MyStructOpt foos[] = {
                { 3.14, 42, "xyz" },
                { 4.49, 5787, "foo" },
                { -3.2, 9873, "zip" },
            };
            MyStructOpt bars[3];

            MPI_Datatype type = mystruct_opt_get_mpi_type();

            cout << "Before MPI_Sendrecv call:\n";
            print_MyStructOpt_arr(bars, 3);

            MPI_Sendrecv(&foos, 3, type, 0, 0, &bars, 3, type, 0, 0, MPI_COMM_SELF, MPI_STATUS_IGNORE);

            cout << "After MPI_Sendrecv call:\n";
            print_MyStructOpt_arr(bars, 3);

            int size; MPI_Type_size(type, &size);
            cout << "Size of MPI_Datatype for MyStruct:  " << size << "\n";
        }

        cout << "Size of MyStruct:  "       << sizeof(MyStruct)     << "\n";
        cout << "Size of MyStructOpt:  "    << sizeof(MyStructOpt)  << "\n";
    }

    MPI_Finalize();
    return 0;
}