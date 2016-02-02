/**
 * @file    main.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the main routine for the CSE6220 programming assignment 1.
 *
 * Copyright (c) 2016 Georgia Institute of Technology. All Rights Reserved.
 */


/*********************************************************************
 *                  !!  DO NOT CHANGE THIS FILE  !!                  *
 *********************************************************************/

#include <mpi.h>
#include <time.h> // for clock_gettime()

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "nqueens.h"
#include "mpi_nqueens.h"

/**
 * Prints the usage of the program.
 */
void print_usage() {
    std::cerr << "Usage: ./mpi_nqueens [options] <n> <k>" << std::endl;
    std::cerr << "      Required arguments:" << std::endl;
    std::cerr << "          <n>     The size of the chess board (i.e. the `n` in n-Queens)." << std::endl;
    std::cerr << "          <k>     The maximum level explored on the master node." << std::endl;
    std::cerr << "      Optional arguments:" << std::endl;
    std::cerr << "          -o      Output all solutions to stdout." << std::endl;
    std::cerr << "          -t      Print tab separated values into one row, the values are" << std::endl;
    std::cerr << "                  (n, k, p, time) in this order." << std::endl;
    std::cerr << "      Example:" << std::endl;
    std::cerr << "          ./mpi_nqueens -o 8 3" << std::endl;
    std::cerr << "                  Will output all solutions to the 8x8 problem where" << std::endl;
    std::cerr << "                  levels [0,k-1] are calculated on the master and [k,n-1]" << std::endl;
    std::cerr << "                  levels are solved on the slave processes." << std::endl;
}


/**
 * @brief Prints all solutions from the local cache.
 */
void print_solutions(const std::vector<unsigned int>& solutions, unsigned int n) {
    size_t num_sols = solutions.size() / n;
    std::cerr << "Printing all " << num_sols << " solutions to stdout:" << std::endl;
    for (size_t i = 0; i < num_sols; ++i) {
        for (unsigned int j = 0; j < n; ++j) {
            if (j != 0)
                std::cout << " ";
            std::cout << solutions[i*n + j];
        }
        std::cout << std::endl;
    }
}

int main(int argc, char *argv[]) {
    // set up MPI
    MPI_Init(&argc, &argv);

    // get communicator size and my rank
    MPI_Comm comm = MPI_COMM_WORLD;
    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    /* code */
    if (rank == 0) {
        // optional arguments
        bool opt_print_solutions = false;
        bool opt_print_table = false;

        // forget about first argument (which is the executable's name)
        argc--;
        argv++;

        // parse optional parameters
        while (argc > 0 && argv[0][0] == '-') {
            char option = argv[0][1];
            switch (option) {
                case 'o':
                    // output all solutions to std out
                    opt_print_solutions = true;
                    break;
                case 't':
                    // print a table row of data
                    opt_print_table = true;
                    break;
                default:
                    print_usage();
                    exit(EXIT_FAILURE);
            }
            // iterate to next argument
            argv++;
            argc--;
        }

        // check that the mandatory parameters are present
        if (argc < 2) {
            print_usage();
            exit(EXIT_FAILURE);
        }
        // parse mandatory parameters
        int n = atoi(argv[0]);
        int k = atoi(argv[1]);
        if (n <= 0 || k <= 0 || k > n) {
            print_usage();
            exit(EXIT_FAILURE);
        }

        // prepare results
        std::vector<unsigned int> results;

        // start timer
        //   we omit the file loading and argument parsing from the runtime
        //   timings, we measure the time needed by the master process
        struct timespec t_start, t_end;
        clock_gettime(CLOCK_MONOTONIC,  &t_start);
        if (p == 1) {
            std::cerr << "[WARNING]: Running the sequential solver. Start with "
                         "mpirun to execute the parallel version." << std::endl;
            // call the sequential solver
            results = nqueens(n);
        } else {
            // call the parallel solver function
            results = master_main(n, k);
        }
        // end timer
        clock_gettime(CLOCK_MONOTONIC,  &t_end);
        // time in seconds
        double time_secs = (t_end.tv_sec - t_start.tv_sec)
                         + (double) (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;

        // print output
        if (opt_print_table) {
            printf("%i\t%i\t%i\t%8.0lf\n", n, k, p, time_secs * 1000.0);
        } else {
            std::cerr << "Number of solutions found: " << results.size()/n << std::endl;
            if (opt_print_solutions) {
                print_solutions(results, n);
            }

            fprintf(stderr, "Run-time of the program: %8.0lf milli-seconds\n", time_secs*1000.0);
        }
    } else {
        worker_main();
    }

    // finalize MPI
    MPI_Finalize();
    return 0;
}
