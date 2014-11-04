/***************************************************************************//**
* \file generateE.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the kernels required to generate elements of the interpolation matrix
*/

#pragma once

#include <types.h>
#include <helpers.h>

/**
* \namespace kernels
* \brief Contain all custom-written CUDA kernels
*/
namespace kernels
{

/**
* \brief Generate the interpolation matrix \c E on the host
*/
void generateEHost(int *ERows,  int *ECols,  real *EVals,
                   int nx, int ny, real *x, real *y, real *dx,
                   int totalPoints, real *xB, real *yB, int *I, int *J);

/**
* \brief Kernel to compute an element of the interpolation matrix \c E on the device
*/
__global__ \
void generateE(int *ERows,  int *ECols,  real *EVals,
               int nx, int ny, real *x, real *y, real *dx,
               int totalPoints, real *xB, real *yB, int *I, int *J);

} // end of namespace kernels
