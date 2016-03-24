/*
    Test program to check that MPI_Datatype for MyStruct works.  To run:

        c++ -std=c++11 single_proc_radix_sort.cpp
        ./a.out
*/

#include <bitset>
#include <numeric>
#include <vector>
#include <iostream>

using namespace std;

#define GET_DIGIT(key, k, offset) ((key) >> (offset)) & ((1 << (k)) - 1)

vector<vector<unsigned int>> create_bucket_table(unsigned int rows, unsigned int columns) {
    return vector<vector<unsigned int>>(rows, vector<unsigned int>(columns));
}

void resize_bucket_table(vector<vector<unsigned int>> &bucket_table, unsigned int rows, unsigned int columns) {
    return bucket_table.resize(rows, vector<unsigned int>(columns));
}

void clean_bucket_table(vector<vector<unsigned int>> &bucket_table) {
    for (auto &row : bucket_table) memset(&row[0], 0, row.size()*sizeof(&row[0]));
    // for (auto &row : bucket_table) std::fill(row.begin(), row.end(), 0);
}

template <typename T>
void prepopulate_bucket_table(vector<vector<unsigned int>> &bucket_table,
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

void print_bucket_table(vector<vector<unsigned int>> &bucket_table) {
    cout << "[ " << bucket_table.size() << " x " <<  bucket_table[0].size() << " ]" << endl;
    for (const auto &row : bucket_table) {
        for (const auto &val : row) {
            cout << val << ' ';
        }  cout << '\n';
    }
}

void perform_prefix_sum_on_bucket_table(vector<vector<unsigned int>> &bucket_table) {
    for (auto &row : bucket_table) std::partial_sum(row.begin(), row.end(), row.begin());
}

template <typename T>
void move_elements_according_to_prefix_sums(vector<T> &lst,
                                            unsigned int (*key_func)(const T&),
                                            vector<vector<unsigned int>> &bucket_table,
                                            unsigned int k,
                                            unsigned int offset) {

    // make copy of array of elements for moving
    auto tmp = lst;

    for (int i=0; i < lst.size(); ++i) {
        unsigned int key = key_func(tmp[i]);
        unsigned int d_i = GET_DIGIT(key, k, offset);

        // cout << lst[i] << "  " << tmp[i] << "  " << key << "  " << d_i << "\n";

        // compute new index of element
        unsigned int new_index = 0;
        for (int j=0; j < d_i; ++j) {
            // cout << bucket_table[j].back() << "+";
            new_index += bucket_table[j].back();

        } // cout << bucket_table[d_i][i-1] << "\n";
        new_index += bucket_table[d_i][i]-1;

        cout << i << " -> " << new_index << "\n";

        // move element over
        lst[new_index] = tmp[i];
    }
}

unsigned int nokey_func(const unsigned int& x) {
    return x;
}

void print_binary(vector<unsigned int> &foo) {
    for (const auto &x : foo) {
        std::bitset<8*sizeof(unsigned int)> bits(x);
        cout << x << ":  " << bits << "\n";
    }
}

template <typename T>
void print_vec(vector<T> &foo) {
    cout << "{ ";
    for (const auto &x : foo) cout << x << ", ";
    cout << "}\n";
}

template <typename T>
void single_processor_radix_sort(vector<T> &lst, unsigned int (*key_func)(const T&), unsigned int k) {
    cout << "Creating bucket table...\n";
    auto bucket_table = create_bucket_table(1 << k, lst.size());
    print_bucket_table(bucket_table);

    unsigned int rnd = 0;
    for (unsigned int offset = 0; offset < 8*sizeof(unsigned int); offset += k, rnd++) {
        cout << "ROUND " << rnd << ": Prepopulating the bucket table...\n";
        prepopulate_bucket_table(bucket_table, lst, nokey_func, k, offset);
        print_bucket_table(bucket_table);

        cout << "ROUND " << rnd << ": Performing prefix sum on bucket table rows...\n";
        perform_prefix_sum_on_bucket_table(bucket_table);
        print_bucket_table(bucket_table);

        move_elements_according_to_prefix_sums(lst, nokey_func, bucket_table, k, offset);
        cout << "LIST AFTER ROUND " << rnd << ":\n";
        print_vec(lst);

        clean_bucket_table(bucket_table);
    }

}

int main(int argc, char *argv[]) {
    vector<unsigned int> lst = {{ 362, 436, 10, 291, 487, 36, 35, 34, 32, 207, 257, 3000, 397 }};
    print_binary(lst);
    single_processor_radix_sort(lst, nokey_func, 3);
    return 0;
}