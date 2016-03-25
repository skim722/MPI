/*
    mpic++ mpi_code_chunks_cpp98.cpp
    mpirun -n 5 a.out
*/

#include <mpi.h>
#include <numeric>
#include <vector>
#include <iostream>

#define GET_DIGIT(key, k, offset) ((key) >> (offset)) & ((1 << (k)) - 1)

using namespace std;

template <typename T>
void print_vec(vector<T> &foo) {
    cout << "{ ";
    for (int i=0; i < foo.size(); ++i) cout << foo[i] << ", ";
    cout << "}\n";
}

vector< vector<unsigned int> > create_bucket_table(unsigned int rows, unsigned int columns) {
    return vector< vector<unsigned int> >(rows, vector<unsigned int>(columns));
}

template <typename T>
void prepopulate_bucket_table(vector< vector<unsigned int> > &bucket_table,
                       vector<T> &lst,
                       unsigned int (*key_func)(const T&),
                       unsigned int k,
                       unsigned int offset) {

    for (int i=0; i < lst.size(); ++i) {
        unsigned int key = key_func(lst[i]);
        unsigned int digit = GET_DIGIT(key, k, offset);
        bucket_table[digit][i] = 1;
    }
}

void perform_prefix_sum_on_bucket_table(vector< vector<unsigned int> > &bucket_table) {
    for (int i=0; i < bucket_table.size(); ++i) {
        vector<unsigned int> &row = bucket_table[i];
        std::partial_sum(row.begin(), row.end(), row.begin());
    }
}


struct tuple_vec {
    vector<unsigned int> A, B, C;

    void print() {
        cout << "-- tuple_vec\n";
        print_vec(A);
        print_vec(B);
        print_vec(C);
        cout << "/- tuple_vec\n";
    }
};

tuple_vec compute_global_t_prime_sums_and_exscans_arrays(const vector< vector<unsigned int> > &prefix_summed_bucket_table,
                                                         MPI_Comm comm) {

    tuple_vec tup;
    vector<unsigned int> &t_primes_summed = tup.A,
                         &t_primes_exscanned = tup.B;

    vector<unsigned int>    t_primes(prefix_summed_bucket_table.size(), 0);
    t_primes_summed         = vector<unsigned int>(prefix_summed_bucket_table.size(), 0);
    t_primes_exscanned      = vector<unsigned int>(prefix_summed_bucket_table.size(), 0);

    // Initialize t_primes with the prefix sums of the bucket table
    for (int i=0; i < t_primes.size(); ++i) {
        t_primes[i] = prefix_summed_bucket_table[i].back();
    }

    // Calculate, for each digit, the total number of elements with the same digit across all processors
    MPI_Allreduce(&t_primes[0], &t_primes_summed[0], t_primes.size(), MPI_UNSIGNED, MPI_SUM, comm);

    // Calculate, for each digit, the total number of elements with the same digit but on processors with smaller rank
    MPI_Exscan(&t_primes[0], &t_primes_exscanned[0], t_primes.size(), MPI_UNSIGNED, MPI_SUM, comm);


    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);
    if (rank == p-1) {
        print_vec(t_primes);
        print_vec(t_primes_summed);
        print_vec(t_primes_exscanned);
    }

    // Return both
    return tup;
}

template <typename T>
tuple_vec compute_G_P_L_EI_arrays(const vector<T> &lst,
                             const vector< vector<unsigned int> > &prefix_summed_bucket_table,
                             const vector<unsigned int> &t_primes_summed,
                             const vector<unsigned int> &t_primes_exscanned,
                             unsigned int (*key_func)(const T&),
                             const unsigned int k,
                             const unsigned int offset) {

    tuple_vec tup;
    vector<unsigned int> &G_EI = tup.A,
                         &P_EI = tup.B,
                         &L_EI = tup.C;

    G_EI = vector<unsigned int>(lst.size(), 0);
    P_EI = vector<unsigned int>(lst.size(), 0);
    L_EI = vector<unsigned int>(lst.size(), 0);

    for (int i=0; i < lst.size(); ++i) {
        // Get the current digit
        unsigned int key = key_func(lst[i]);
        unsigned int d_i = GET_DIGIT(key, k, offset);

        // Sum up the T''s to get the total number of elements with digit smaller than the current digit
        for (int j=0; j < d_i; ++j) {
            G_EI[i] += t_primes_summed[j];
        }

        // Get the number of elements with the same digit, but on processors with smaller rank
        P_EI[i] = t_primes_exscanned[d_i];

        // Get the number of elements with the same digit on the same processor and LEFT of this element (hence the -1)
        prefix_summed_bucket_table[d_i][i] - 1;
    }

    return tup;
}

template <typename T>
vector<unsigned int> compute_global_target_index(const vector<T> &lst,
                                                 const vector< vector<unsigned int> > &prefix_summed_bucket_table,
                                                 unsigned int (*key_func)(const T&),
                                                 MPI_Comm comm,
                                                 unsigned int k,
                                                 unsigned int offset) {

    // Compute the T'' sums and exscans (requires MPI)
    tuple_vec global_t_prime_sums_and_exscans = compute_global_t_prime_sums_and_exscans_arrays(prefix_summed_bucket_table, comm);

    // Using the T' sums and exscans, compute G_EI, P_EI, and L_EI
    tuple_vec index_arrays = compute_G_P_L_EI_arrays(lst,
                                                     prefix_summed_bucket_table,
                                                     global_t_prime_sums_and_exscans.A,
                                                     global_t_prime_sums_and_exscans.B,
                                                     key_func, k, offset);

    // Sum them together to get the full T_EI
    vector<unsigned int> T_EI(index_arrays.A.size(), 0);
    for (int i=0; i < index_arrays.A.size(); ++i) {
        T_EI[i] = index_arrays.A[i] + index_arrays.B[i] + index_arrays.C[i];
    }

    return T_EI;
}

unsigned int nokey_func(const unsigned int& x) {
    return x;
}

void print_bucket_table(vector< vector<unsigned int> > &bucket_table) {
    cout << "[ " << bucket_table.size() << " x " <<  bucket_table[0].size() << " ]" << endl;
    for (int i=0; i < bucket_table.size(); ++i) {
        vector<unsigned int> &row = bucket_table[i];
        for (int j=0; j < row.size(); ++j) {
            cout << row[j] << ' ';
        }  cout << '\n';
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);


    vector<unsigned int> lst(10, 2);
    unsigned int k = 3, offset=0;

    vector< vector<unsigned int> > bucket_table = create_bucket_table(1 << k, lst.size());
    prepopulate_bucket_table(bucket_table, lst, nokey_func, k, offset);
    perform_prefix_sum_on_bucket_table(bucket_table);

    if (rank == 0) {
        print_bucket_table(bucket_table);
    }

    vector<unsigned int> global_indexes = compute_global_target_index(lst, bucket_table, nokey_func, comm, k, offset);

    if (rank == p-1) {
        print_vec(global_indexes);
    }

    MPI_Finalize();

    return 0;
}
