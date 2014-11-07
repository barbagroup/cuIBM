/***************************************************************************//**
 * \file generateQT.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the kernels to generate gradient matrix and interpolation matrix.
 */


#pragma once

#include <types.h>
#include <helpers.h>


/**
 * \namespace kernels
 * \brief Contains all custom-written CUDA kernels.
 */
namespace kernels
{

__global__ \
void updateQ(int *QRows, int *QCols, real *QVals, int QSize, int *tags);

// generate the divergence matrix
void generateQT(int *QTRows, int *QTCols, real *QTVals, int nx, int ny);

// update elements of the divergence matrix and the interpolation matrix
__global__ \
void updateQT(int *QTRows, int *QTCols, real *QTVals,
              int *ERows,  int *ECols,  real *EVals,
              int nx, int ny, real *x, real *y, real *dx,
              int totalPoints, real *xB, real *yB, int *I, int *J);

// update the divergence matrix and the interpolation matrix
void updateQTHost(int *QTRows, int *QTCols, real *QTVals,
              int *ERows,  int *ECols,  real *EVals,
              int nx, int ny, real *x, real *y, real *dx,
              int totalPoints, real *xB, real *yB, int *I, int *J);

} // end of namespace kernels
