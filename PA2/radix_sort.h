/**
 * @file    radix_sort.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 *
 * Copyright (c) 2016 Georgia Institute of Technology. All Rights Reserved.
 */

/*
 * TODO: implement your radix sort solution in this file
 */

#include <mpi.h>
#include <numeric>
#include <vector>
#include <iostream>

// returns the value of the digit starting at offset `offset` and containing `k` bits
#define GET_DIGIT(key, k, offset) ((key) >> (offset)) & ((1 << (k)) - 1)

using namespace std;

int get_rank(MPI_Comm comm) {
    int rank; MPI_Comm_rank(comm, &rank); return rank;
}

template <typename T>
void print_vec(vector<T> &foo, unsigned int (*key_func)(const T&)) {
    cout << "{ ";
    for (int i=0; i < foo.size(); ++i) cout << key_func(foo[i]) << ", ";
    cout << "}\n";
}

template <typename T>
void mpi_print_vec(vector<T> &lst, unsigned int (*key_func)(const T&), MPI_Comm comm) {
    int num_processors, rank;
    MPI_Comm_size(comm, &num_processors);
    MPI_Comm_rank(comm, &rank);

    MPI_Barrier(comm);
    for (int i=0; i < num_processors; ++i) {
        if (i == rank) {
            cout << "PROCESSOR [" << i << "]:  ";
            print_vec(lst, key_func);
            cout << "\n";
        }
        MPI_Barrier(comm);
    }
}

struct tuple_vec {
    vector<unsigned int> A, B, C;
};

vector< vector<unsigned int> > create_bucket_table(unsigned int rows, unsigned int columns) {
    return vector< vector<unsigned int> >(rows, vector<unsigned int>(columns));
}

void clean_bucket_table(vector< vector<unsigned int> > &bucket_table) {
    for (int i=0; i < bucket_table.size(); ++i) {
        std::fill(bucket_table[i].begin(), bucket_table[i].end(), 0);
    }
}

void print_bucket_table(const vector< vector<unsigned int> > &bucket_table) {
    cout << "[ " << bucket_table.size() << " x " <<  bucket_table[0].size() << " ]" << endl;
    for (int i=0; i < bucket_table.size(); ++i) {
        const vector<unsigned int> &row = bucket_table[i];
        for (int j=0; j < row.size(); ++j) {
            cout << row[j] << ' ';
        }  cout << '\n';
    }
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

template <typename T>
void move_elements_according_to_prefix_sums(vector<T> &lst,
                                            unsigned int (*key_func)(const T&),
                                            vector< vector<unsigned int> > &prefix_summed_bucket_table,
                                            unsigned int k,
                                            unsigned int offset) {
    // make a temporary copy of array of elements for moving
    vector<T> tmp = lst;
    for (int i=0; i < lst.size(); ++i) {
        unsigned int key = key_func(tmp[i]);
        unsigned int d_i = GET_DIGIT(key, k, offset);

        // compute new index of element
        unsigned int new_index = 0;
        for (int j=0; j < d_i; ++j) {
            new_index += prefix_summed_bucket_table[j].back();
        }

        new_index += prefix_summed_bucket_table[d_i][i]-1;

        // move element over
        // cout << i << " -> " << new_index << "\n";
        lst[new_index] = tmp[i];
    }
}

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
        L_EI[i] = prefix_summed_bucket_table[d_i][i] - 1;
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

template <typename T>
void exchange_elements_between_processors(vector<T> &lst,
                                          const vector<unsigned int> &global_indexes,
                                          MPI_Datatype mpi_dt,
                                          MPI_Comm comm) {
    // Get number of processors
    int num_processors; MPI_Comm_size(comm, &num_processors);

    // Run histogram over the list of elements to compute how many items will be sent to how which processors
    vector<int> send_counts(num_processors, 0);
    for (int i=0; i < global_indexes.size(); ++i) {
        // global_indexes.size() SHOULD BE EQUAL to lst.size() (it is a BUG if it's not!)
        int destined_processor = global_indexes[i] / lst.size();
        send_counts[destined_processor]++;
    }

    // Initialize receive_counts
    vector<int> receive_counts(num_processors, 0);

    // Exchange the sendcount information with the other processors so a receive_counts array can be filled in for MPI_Alltoallv use later on
    // (Tell the other processors how much data will be coming)
    MPI_Alltoall(&send_counts[0], 1, MPI_INT, &receive_counts[0], 1, MPI_INT, comm);

    // Initialize displacement arrays for the send and receive buffers
    vector<int> send_displacements(num_processors, 0), receive_displacements(num_processors, 0);
    for (int i=1; i < num_processors; ++i) {
        // calculate displacements
        send_displacements[i]       = send_displacements[i-1] + send_counts[i-1];
        receive_displacements[i]    = receive_displacements[i-1] + receive_counts[i-1];
    }

    // calculate receive size (should be == lst.size())
    int receive_size = std::accumulate(receive_counts.begin(), receive_counts.end(), 0);

    // Initialize the receive buffer
    vector<T> recvbuf(receive_size);

    // send/receive different amounts of data to/from each processor
    MPI_Alltoallv(&lst[0], &send_counts[0], &send_displacements[0], mpi_dt,
                  &recvbuf[0], &receive_counts[0], &receive_displacements[0], mpi_dt, comm);

    // Copy revbuf over to lst
    lst = recvbuf;
    return;
}










/**
 * @brief   Parallel distributed radix sort.
 *
 * This function sorts the distributed input range [begin, end)
 * via lowest-significant-byte-first radix sort.
 *
 * This function will sort elements of type `T`, by the key of type `unsigned int`
 * which is returned by the key-access function `key_func`.
 *
 * The MPI datatype for the templated (general) type `T` has to be passed
 * via the `dt` parameter.
 *
 * @param begin         A pointer to the first element in the local range to be sorted.
 * @param end           A pointer to the end of the range to be sorted. This
 *                      pointer points one past the last element, such that the
 *                      total number of elements is given by `end - begin`.
 * @param key_func      A function with signature: `unsigned int (const T&)`.
 *                      This function returns the key of each element, which is
 *                      used for sorting.
 * @param dt            The MPI_Datatype which represents the type `T`. This
 *                      is used whenever elements of type `T` are communicated
 *                      via MPI.
 * @param comm          The communicator on which the sorting happens.
 *                      NOTE: this is not necessarily MPI_COMM_WORLD. Call
 *                            all MPI functions with this communicator and
 *                            NOT with MPI_COMM_WORLD.
 */
template <typename T>
void radix_sort(T* begin, T* end, unsigned int (*key_func)(const T&), MPI_Datatype dt, MPI_Comm comm, unsigned int k = 16) {
    // Copy the input array into vector for easier manipulation
    vector<T> lst(begin, end);

    // the number of histogram buckets = 2^k = 1 << k
    vector< vector<unsigned int> > bucket_table = create_bucket_table(1 << k, lst.size());

    for (unsigned int offset = 0; offset < 8*sizeof(unsigned int); offset += k) {
        // 1.) create histogram and sort via bucketing (~ counting sort)
        clean_bucket_table(bucket_table);
        prepopulate_bucket_table(bucket_table, lst, key_func, k, offset);
        perform_prefix_sum_on_bucket_table(bucket_table);
        move_elements_according_to_prefix_sums(lst, key_func, bucket_table, k, offset);


        // 2.) get global histograms (P, G) via MPI_Exscan/MPI_Allreduce,...
                // We must first reset the bucket table before proceeding, because the elements have been moved around
        clean_bucket_table(bucket_table);
        prepopulate_bucket_table(bucket_table, lst, key_func, k, offset);
        perform_prefix_sum_on_bucket_table(bucket_table);
        vector<unsigned int> global_indexes = compute_global_target_index(lst, bucket_table, key_func, comm, k, offset);


        // 3.) calculate send_counts
        // 4.) communicate send_counts to get recv_counts
        // 4.) calculate displacements
        // 6.) MPI_Alltoallv
        exchange_elements_between_processors(lst, global_indexes, dt, comm);


        // 7.) local sorting via bucketing (~ counting sort)
        clean_bucket_table(bucket_table);
        prepopulate_bucket_table(bucket_table, lst, key_func, k, offset);
        perform_prefix_sum_on_bucket_table(bucket_table);
        move_elements_according_to_prefix_sums(lst, key_func, bucket_table, k, offset);
    }

    // Copy results back to the original array that was passed in
    std::copy(lst.begin(), lst.end(), begin);
}

