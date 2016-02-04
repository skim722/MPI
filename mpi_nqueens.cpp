/**
 * @file    mpi_nqueens.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the parallel, master-worker nqueens solver.
 *
 * Copyright (c) 2016 Georgia Institute of Technology. All Rights Reserved.
 */

/*********************************************************************
 *                  Implement your solutions here!                   *
 *********************************************************************/

#include "mpi_nqueens.h"

#include <mpi.h>
#include <vector>
#include "nqueens.h"

#define TERMINATION_MSG 42

/**
 * @brief The master's call back function for each found solution.
 *
 * This is the callback function for the master process, that is called
 * from within the nqueens solver, whenever a valid solution of level
 * `k` is found.
 *
 * This function will send the partial solution to a worker which has
 * completed his previously assigned work. As such this function must
 * also first receive the solution from the worker before sending out
 * the new work.
 *
 * @param solution      The valid solution. This is passed from within the
 *                      nqueens solver function.
 */
void master_solution_func(std::vector<unsigned int>& solution) {
    // TODO: receive solutions or work-requests from a worker and then
    //       proceed to send this partial solution to that worker.

    // Needs to be a NON-blocking receive/send
    // MPI_Recv(void* buf, int count, MPI_Datatype type, int source, int tag, MPI_Comm comm, MPI_Status* status);
    // MPI_Send(void* buf, int count, MPI_Datatype type, int dest, int tag, MPI_Comm comm);

    // http://stackoverflow.com/questions/10490983/mpi-slave-processes-hang-when-there-is-no-more-work
}



/**
 * @brief   Performs the master's main work.
 *
 * This function performs the master process' work. It will sets up the data
 * structure to save the solution and call the nqueens solver by passing
 * the master's callback function.
 * After all work has been dispatched, this function will send the termination
 * message to all worker processes, receive any remaining results, and then return.
 *
 * @param n     The size of the nqueens problem.
 * @param k     The number of levels the master process will solve before
 *              passing further work to a worker process.
 */
std::vector<unsigned int> master_main(unsigned int n, unsigned int k) {
    // TODO: send parameters (n,k) to workers via broadcast (MPI_Bcast)
    // MPI_Bcast(void* data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator)
    MPI_Bcast(&n, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    // allocate the vector for the solution permutations
    std::vector<unsigned int> pos(n);

    // generate all partial solutions (up to level k) and call the
    // master solution function
    nqueens_by_level(pos, 0, k, &master_solution_func);

    // TODO: get remaining solutions from workers and send termination messages

    // TODO: return all combined solutions
}


// stores all local solutions.
struct WorkerSolutionStore {
    // store solutions in a static member variable
    static std::vector<unsigned int>& solutions() {
        static std::vector<unsigned int> sols;
        return sols;
    }
    static void add_solution(const std::vector<unsigned int>& sol) {
        // add solution to static member
        solutions().insert(solutions().end(), sol.begin(), sol.end());
    }
    static void clear_solutions() {
        solutions().clear();
    }
};

/**
 * @brief The workers' call back function for each found solution.
 *
 * This is the callback function for the worker processes, that is called
 * from within the nqueens solver, whenever a valid solution is found.
 *
 * This function saves the solution into the worker's solution cache.
 *
 * @param solution      The valid solution. This is passed from within the
 *                      nqueens solver function.
 */
void worker_solution_func(std::vector<unsigned int>& solution) {
    // TODO: save the solution into a local cache
    WorkerSolutionStore::add_solution(solution);
}

/**
 * @brief   Performs the worker's main work.
 *
 * This function implements the functionality of the worker process.
 * The worker will receive partially completed work items from the
 * master process and will then complete the assigned work and send
 * back the results. Then again the worker will receive more work from the
 * master process.
 * If no more work is available (termination message is received instead of
 * new work), then this function will return.
 */
void worker_main() {
    unsigned int n, k;

    // Get rank
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // TODO receive the parameters `n` and `k` from the master process via MPI_Bcast
    MPI_Bcast(&n, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    // TODO: implement the worker's functions: receive partially completed solutions,
    //       calculate all possible solutions starting with these queen positions
    //       and send solutions to the master process. then ask for more work.
    while (true) {
        // Prepare the buffer
        std::vector<unsigned int> partial_solution(n);
        MPI_Status stat;

        // Receive the buffer of size n and begin computation
        // MPI_Recv(void* data, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm communicator, MPI_Status* status)
        MPI_Recv(&partial_solution[0], n, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

        // If we get the termination signal, exit the loop
        if (stat.MPI_TAG == TERMINATION_MSG) {
            printf("Process %d got the termination signal; exiting work loop.\n", rank);
            break;
        }

        // Else, compute (saves solution to a long vector cache)
        nqueens_by_level(partial_solution, k, n, &worker_solution_func);
        std::vector<unsigned int> &vec_buffer = WorkerSolutionStore::solutions();
        unsigned int buffer_size = vec_buffer.size();

        // Send the solutions buffer size (m solutions means the solutions vector is of size m x n), so the recipient knows how much to receive
        MPI_Send(&buffer_size, 1, MPI_UNSIGNED, 0, 528492, MPI_COMM_WORLD);
        MPI_Send(&vec_buffer[0], buffer_size, MPI_UNSIGNED, 0, 528493, MPI_COMM_WORLD);
    }
}


