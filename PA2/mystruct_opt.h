/**
 * @file    mystruct_opt.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 *
 * Copyright (c) 2016 Georgia Institute of Technology. All Rights Reserved.
 */
#ifndef MYSTRUCT_H
#define MYSTRUCT_H

#include <mpi.h>
#include <iostream>


/*
 * TODO: Investigate how you can change this structure to increase MPI
 *       performance within radix sort.
 */
struct MyStruct {
    // Order the struct member variables as such to remove 8 bytes of local padding (reduce struct offset (sizeof()) from 24 to 16 bytes)
    double d;           // 8 bytes
    unsigned int key;   // 4 bytes
    char e[4];          // 4 bytes
};

/**
 * Returns the MPI_Datatype for `MyStruct`.
 *
 * TODO: You have to implement this function in mystruct_opt.cpp
 */
MPI_Datatype mystruct_get_mpi_type();


/**
 * Returns a random instance of `MyStruct`.
 */
MyStruct mystruct_rand();

/**
 * Returns the value of the key of `MyStruct`.
 */
unsigned int mystruct_key_access(const MyStruct& s);


/**
 * Further utility functions that you don't have to worry about.
 */

// operator used for comparison based sorting
bool operator<(const MyStruct& lhs, const MyStruct& rhs);
// operator for comparing if two MyStruct's are equal
bool operator==(const MyStruct& lhs, const MyStruct& rhs);
// output format for MyStruct
std::ostream& operator <<(std::ostream& stream, const MyStruct& s);


#endif
