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
      int rank;
      MPI_Init(&argc,&argv);
      MPI_Comm_rank (MPI_COMM_WORLD,&rank);
      if(rank==0){
         MPI_BCast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
      }
    // allocate the vector for the solution permutations
    std::vector<unsigned int> pos(n);

    // generate all partial solutions (up to level k) and call the
    // master solution function
    nqueens_by_level(pos, 0, k, &master_solution_func);

    // TODO: get remaining solutions from workers and send termination messages

    // TODO: return all combined solutions
    return MPI_Finalize();
}

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
    // TODO receive the parameters `n` and `k` from the master process via MPI_Bcast

    // TODO: implement the worker's functions: receive partially completed solutions,
    //       calculate all possible solutions starting with these queen positions
    //       and send solutions to the master process. then ask for more work.
}


