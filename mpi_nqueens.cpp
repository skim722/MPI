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

#define WORKER_READY            40
#define NEW_WORK_MSG            41
#define TERMINATION_MSG         42
#define WORK_DONE_MSG_N         43
#define WORK_DONE_MSG_BUFFER    44


// stores all local solutions for the worker
struct MasterSolutionStore {
    // store solutions in a static member variable
    static std::vector<unsigned int>& solutions() {
        static std::vector<unsigned int> sols;
        return sols;
    }
    static void add_solutions(const std::vector<unsigned int>& sol) {
        // add solution to static member
        solutions().insert(solutions().end(), sol.begin(), sol.end());
    }
    static void clear_solutions() {
        solutions().clear();
    }
};

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
    unsigned int buffer_size;
    MPI_Status stat;

    // SK
    // if (receive sol(1,n) from proc p) Send partial_sol(1,k) to p
    // else if(receive work request from proc p) Send partial_sol(1,k) to p

    // First receive a message from any source.
    // By the symmetry of the code, the first message will either be WORKER_READY or WORK_DONE_MSG_N, and buffer will be size 1
    MPI_Recv(&buffer_size, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    printf("Process 0 got tag %d\n", stat.MPI_TAG);

    // If the tag is WORK_DONE_MSG_N, then we must first receive the finished results before sending new work over
    // Otherwise, we continue (buffer_size won't be needed)
    if (stat.MPI_TAG == WORK_DONE_MSG_N) {
        std::vector<unsigned int> finished_solutions(buffer_size);

        // Receive the actual solution set, from the SAME worker (stat.MPI_SOURCE)
        MPI_Recv(&finished_solutions[0], buffer_size, MPI_UNSIGNED, stat.MPI_SOURCE, WORK_DONE_MSG_BUFFER, MPI_COMM_WORLD, &stat);
        printf("Process 0 got tag %d from process %d\n", stat.MPI_TAG, stat.MPI_SOURCE);
printf("buffer size=%d\n", buffer_size);
        //printf("                  !!sol 0 to 6 : %d %d %d | %d %d \n",finished_solutions[0],finished_solutions[1],finished_solutions[2],finished_solutions[3],finished_solutions[4],finished_solutions[5]);
        // Store the completed results that came from the worker
        MasterSolutionStore::add_solutions(finished_solutions);
    
    }
    // Send new work to the SAME worker
    printf("Process 0 sends tag %d to process %d\n", NEW_WORK_MSG, stat.MPI_SOURCE);
  //  printf("                  sol 0 to 2 : %d %d %d | %d %d \n",solution[0],solution[1],solution[2],solution[3],solution[4],solution[5]);
    MPI_Send(&solution[0], solution.size(), MPI_UNSIGNED, stat.MPI_SOURCE, NEW_WORK_MSG, MPI_COMM_WORLD);
    //}
    // SK
    // http://stackoverflow.com/questions/10490983/mpi-slave-processes-hang-when-there-is-no-more-work

    /*MPI_Status stat;
    for(int i = 1; i < num_procs; i++){
        MPI_Recv(&buffer_size, 1, MPI_UNSIGNED, i, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
        if (stat.MPI_TAG == 528492) {
             MPI_Send(&partial_solution[0], n, MPI_UNSIGNED, i, MPI_ANY_TAG, MPI_COMM_WORLD);

        }
    }*/
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
    int num_processors;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
    for (int w = 1; w < num_processors; ++w) {
        MPI_Status stat;
        unsigned int buffer_size;
        MPI_Recv(&buffer_size, 1, MPI_UNSIGNED, w, WORK_DONE_MSG_N, MPI_COMM_WORLD, &stat);

        // Receive the actual solution set, from the SAME worker (stat.MPI_SOURCE)
        std::vector<unsigned int> finished_solutions(buffer_size);
        MPI_Recv(&finished_solutions[0], buffer_size, MPI_UNSIGNED, w, WORK_DONE_MSG_BUFFER, MPI_COMM_WORLD, &stat);

        // Store the completed results that came from the worker
        MasterSolutionStore::add_solutions(finished_solutions);

        // Send termination message
        MPI_Send(&pos[0], pos.size(), MPI_UNSIGNED, w, TERMINATION_MSG, MPI_COMM_WORLD);
    }

    // TODO: return all combined solutions
    return MasterSolutionStore::solutions();
}


// stores all local solutions for the worker
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

    // Send the first work request message (the tag is the important value to send over, n is just a filler value)
    MPI_Send(&n, 1, MPI_UNSIGNED, 0, WORKER_READY, MPI_COMM_WORLD);

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
        printf("Process slave got tag %d\n", stat.MPI_TAG);

        // If we get the termination signal, exit the loop
        // http://stackoverflow.com/questions/10490983/mpi-slave-processes-hang-when-there-is-no-more-work
        if (stat.MPI_TAG == TERMINATION_MSG) {
            printf("Process %d got the termination signal; exiting work loop.\n", rank);
            break;
        }

        // Else, compute (saves solution to a long vector cache)
        nqueens_by_level(partial_solution, k, n, &worker_solution_func);
        std::vector<unsigned int> &vec_buffer = WorkerSolutionStore::solutions();
        unsigned int buffer_size = vec_buffer.size();

        // Send the solutions buffer size (m solutions means the solutions vector is of size m x n), so the recipient knows how much to receive
        printf("Process %d sends tag %d\n", rank, WORK_DONE_MSG_N);
        MPI_Send(&buffer_size, 1, MPI_UNSIGNED, 0, WORK_DONE_MSG_N, MPI_COMM_WORLD);
        printf("Process %d sends tag %d\n", rank, WORK_DONE_MSG_BUFFER);
        MPI_Send(&vec_buffer[0], buffer_size, MPI_UNSIGNED, 0, WORK_DONE_MSG_BUFFER, MPI_COMM_WORLD);
        // MPI_Send(&buffer_size, 1, MPI_UNSIGNED, 0, 528492, MPI_COMM_WORLD);
        // MPI_Send(&vec_buffer[0], buffer_size, MPI_UNSIGNED, 0, 528493, MPI_COMM_WORLD); //ask more work tag: 528492 528493
    }
}


