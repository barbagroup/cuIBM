/***************************************************************************//**
* \file generateA.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Contain declaration of the kernels required to generate the matrix A
*/

#pragma once

#include <types.h>

/**
* \namespace kernels
* \brief Contain all the custom-written CUDA kernels
*/
namespace kernels
{

/**
* \brief Generate the matrix A
*/
__global__
void generateA(int *ARows, int *ACols, real *AVals, real *MVals, int *LRows, int *LCols, real *LVals, int ASize, real alpha);

/**
* \brief Generate the matrix A for the direct forcing method
*/
__global__
void generateADirectForcing(int *ARows, int *ACols, real *AVals, real *MVals, int *LRows, int *LCols, real *LVals, int ASize, real alpha, int *tags);

} // end of namespace kernels
