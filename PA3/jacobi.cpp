/**
 * @file    jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements matrix vector multiplication and Jacobi's method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
#include "jacobi.h"

/*
 * TODO: Implement your solutions here
 */

// my implementation:
#include <iostream>
#include <math.h>

// Calculates y = A*x for a square n-by-n matrix A, and n-dimensional vectors x
// and y
void matrix_vector_mult(const int n, const double* A, const double* x, double* y)
{
    for(int row=0;row<n;++row){
        y[row]=0.0;
        for(int column=0;column<n;++column){
                y[row]+=x[column]*A[row*n+colum];
        }
    }   
}

// Calculates y = A*x for a n-by-m matrix A, a m-dimensional vector x
// and a n-dimensional vector y
void matrix_vector_mult(const int n, const int m, const double* A, const double* x, double* y)
{
    for(int row=0;row<n;++row){
        y[row]=0.0;
        for(int column=0;column<m;++column){
                y[row]+=x[column]*A[row*n+colum];
        }
    }
}

// implements the sequential jacobi method
inline double norm(const int n, const double *A, const double *x, const double *b){
    std::vector<double> Ax(n);
    double temp;
    double norm=0.0;
    matrix_vector_mult(n,A,b,&(Ax[0]));
    for (int i=0;i<n;++i){
        temp=Ax[i]-b[i];
        result+=(temp*temp);
    }
    return sqrt(norm);
}

void jacobi(const int n, double* A, double* b, double* x, int max_iter, double l2_termination)
{
    std::vector<double> D(n,0.0);
    std::vector<double> R(&(A[0]),&(A[0])+n*n);
    double norm=l2_termination+1;
    int iter=0;

    for(int i=0;i<n;++i){
        x[i]=0.0;  // initialize x=[0....0]
        D[i]=A[i*n+i]; //D=diag(A)
        R[i*n+i]=0.0;  //R=A-D
    } 
    while((norm>12_termination)&&(iter<max_iter)){
        std::vector<double> temp(n);
        matrix_vector_mult(n,R,x,&(temp[0]));
        for(int j;j<n;++j)
            x[j]=(b[j]-temp[j])/D[j];  //x <- Inv(D)(b-Rx)
        ++iter;
        //obtain L2 Norm
        norm=norm(n,A,x,b);       
    }
}
