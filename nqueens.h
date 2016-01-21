/**
 * @file    nqueens.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Declares the functions for solving the nqueens problem in.
 *
 * Copyright (c) 2016 Georgia Institute of Technology. All Rights Reserved.
 */

/*********************************************************************
 *                  !!  DO NOT CHANGE THIS FILE  !!                  *
 *********************************************************************/

#ifndef NQUEENS_H
#define NQUEENS_H

#include <vector>

/**
 * @brief   Generates the solutions for the n-queen problem in a specified range.
 *
 *
 * This function will search for all valid solutions in the range
 * [start_level,max_level). This assumes that the partial solution for levels
 * [0,start_level) is a valid solution.
 *
 * The master and the workers will both use this function with different parameters.
 * The master will call the function with start_level=0 and max_level=k, and the
 * workers will call this function with start_level=k and max_level=n.
 *
 * for every valid solution that is found, the given callback function is called
 * with the current valid solution as parameter. Depending on what is returned
 * by the callback function, this function will either quit or continue to
 * generate more solutions. This return value is used to either generate
 * all possible solutions or quit after the first valid solution has been
 * found on one worker.
 *
 * @param pos   A vector of length n that contains the partial solution for the
 *              range [0, start_level).
 * @param start_level   Gives the level at which the algorithm will start to
 *                      generate solutions. The solution vector `pos` is assumed
 *                      to have a valid solution for positions [0,start_level-1].
 * @param max_level     The level at which the algorithm will stop and pass
 *                      the current valid solution to the callback function.
 * @param success_func  A callback function that is called whenever a valid
 *                      solution of `max_level` levels is found. The found
 *                      solution is passed as parameter to this callback function.
 */
void nqueens_by_level(std::vector<unsigned int> pos, unsigned int start_level,
                      unsigned int max_level,
                      void (* const success_func)(std::vector<unsigned int>&));

/**
 * @brief   Returns all solutions for the n-queens problem.
 *
 * This will calculate and return all solutions to the problem.
 * Each solution is represented as `n` integers. Solutions are concatenated
 * together, such that the returned vector will contain m*n integers,
 * if there are `m` solutions in total.
 *
 * @param n     The `n` in n-queens. This is the size of the chessboard for
 *              which to solve the problem.
 * @returns     All solutions to the n-queens problem. 
 */
std::vector<unsigned int> nqueens(unsigned int n);

#endif // NQUEENS_H
