


std::tuple<vector<unsigned int>, vector<unsigned int>> compute_global_t_prime_sums_and_exscans_arrays(const vector<vector<unsigned int>> &prefix_summed_bucket_table,
                                                                                                      MPI_Comm comm) {
    vector<unsigned int> t_primes(prefix_summed_bucket_table.size()),
                         t_primes_summed(prefix_summed_bucket_table.size()),
                         t_primes_exscanned(prefix_summed_bucket_table.size());

    // Initialize t_primes with the prefix sums of the bucket table
    for (auto i=0; i < t_primes.size(); ++i) {
        t_primes[i] = prefix_summed_bucket_table[i].back();
    }

    // Calculate, for each digit, the total number of elements with the same digit across all processors
    MPI_Allreduce(&t_primes[0], &t_primes_summed[0], t_primes.size(), MPI_UNSIGNED, MPI_SUM, comm);

    // Calculate, for each digit, the total number of elements with the same digit but on processors with smaller rank
    MPI_Exscan(&t_primes[0], &t_primes_exscanned[0], t_primes.size(), MPI_UNSIGNED, MPI_SUM, comm);

    // Return both
    return make_tuple(t_primes_summed, t_primes_exscanned);
}

template <typename T>
std::tuple<vector<unsigned int>, vector<unsigned int>, vector<unsigned int>> compute_G_P_L_EI_arrays(const vector<T> &lst,
                                                                                                     const vector<vector<unsigned int>> &prefix_summed_bucket_table,
                                                                                                     const vector<unsigned int> &t_primes_summed,
                                                                                                     const vector<unsigned int> &t_primes_exscanned,
                                                                                                     unsigned int (*key_func)(const T&),
                                                                                                     const unsigned int k,
                                                                                                     const unsigned int offset) {

    vector<unsigned int> G_EI(lst.size(), 0),
                         P_EI(lst.size(), 0),
                         L_EI(lst.size(), 0);

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

    return make_tuple(G_EI, P_EI, L_EI);
}

vector<unsigned int> compute_global_target_index(const vector<T> &lst,
                                                 const vector<vector<unsigned int>> &prefix_summed_bucket_table,
                                                 unsigned int (*key_func)(const T&),
                                                 MPI_Comm comm,
                                                 unsigned int k,
                                                 unsigned int offset) {

    // Compute the T'' sums and exscans (requires MPI)
    vector<unsigned int> t_primes_summed, t_primes_exscanned;
    std::tie(t_primes_summed, t_primes_exscanned) = compute_global_t_prime_sums_and_exscans_arrays(prefix_summed_bucket_table, comm);

    // Using the T' sums and exscans, compute G_EI, P_EI, and L_EI
    vector<unsigned int> G_EI, P_EI, L_EI;
    std::tie(G_EI, P_EI, L_EI) = compute_G_P_L_EI_arrays(lst, prefix_summed_bucket_table, t_primes_summed, t_primes_exscanned, key_func, k, offset);

    // Sum them together to get the full T_EI
    vector<unsigned int> T_EI(G_EI.size());
    for (int i=0; i < G_EI.size(); ++i) {
        T_EI[i] = G_EI[i] + P_EI[i] + L_EI[i];
    }

    return T_EI;
}