/**
 * @file    mpi_nqueens.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Declares the functions for solving the nqueens problem in
 *          parallel using a master-worker paradigm.
 *
 * Copyright (c) 2016 Georgia Institute of Technology. All Rights Reserved.
 */

/*********************************************************************
 *                  !!  DO NOT CHANGE THIS FILE  !!                  *
 *********************************************************************/

#ifndef MPI_NQUEENS_H
#define MPI_NQUEENS_H

#include <vector>

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
std::vector<unsigned int> master_main(unsigned int n, unsigned int k);

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
void worker_main();

#endif // MPI_NQUEENS_H
