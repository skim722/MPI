CSE 6220 Programming Assignment 1
=================================

We provide this framework to get you started. There is a `main.cpp` file, which
is the main entry point to the program you have to write. This is already
fully implemented. This function reads command line parameters, reads the
input file and finally writes the output file. You should take a look at this
file, but you won't have to change anything in it.

For this programming assignment, you will have to implement the functions
that are in the files `nqueens.cpp`, and `mpi_nqueens.cpp`.
The functions that you have to implement are declared and documented in the
corresponding `.h` files.

Compiling
---------
We provide a makefile for compiling the code. Simply run `make` in the
projects directory to build the executable: `nqueens` and
Run this executable without command-line parameters for seeing usage
descriptions.

Running on Jinx
---------------
To run your parallel program on jinx, we supply a very basic PBS script
in `pbs_script.sh`. You'll have to change this script to point to the correct
project folder. You'll also have to modify the script to achieve multiple runs
for your experiments and for testing different inputs. DO NOT just take this
script as it is.

Sample Output
-------------
We give you some sample outputs in the `sample` folder (for `n=8, n=10, n=13`),
which you can use as a reference to compare if your solutions are correct. For
grading, we will compare the output of your program with these and other, much
larger instances.

The order of solutions can be arbitrary, so for comparison, you can use the
linux tool `sort` and then `diff`:

```sh
    diff <(sort your_8.txt) <(sample/8.txt)
```
