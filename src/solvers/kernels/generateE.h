/**
 * \file generateE.h
 * \brief Declaration of the kernels to generate the interpolation matrix.
 */


#pragma once

#include "utilities/types.h"
#include "utilities/helpers.h"


/**
 * \namespace kernels
 * \brief Contains all custom-written CUDA kernels.
 */
namespace kernels
{

// generate the interpolation matrix (on the host)
void generateEHost(int *ERows,  int *ECols,  real *EVals,
                   int nx, int ny, real *x, real *y, real *dx,
                   int totalPoints, real *xB, real *yB, int *I, int *J);

// compute elements of the interpolation matrix
__global__ \
void generateE(int *ERows,  int *ECols,  real *EVals,
               int nx, int ny, real *x, real *y, real *dx,
               int totalPoints, real *xB, real *yB, int *I, int *J);

} // end of namespace kernels
