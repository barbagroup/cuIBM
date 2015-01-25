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
// update the Q matrix for the direct forcing method
__global__ \
void updateQ(int *QRows, int *QCols, real *QVals, int QSize, int *tags);

// update the QT matrix for the direct forcing method
__global__ \
void updateQT(int *QTRows, int *QTCols, real *QTVals, int QTSize, int *tags, real *coeffs);

// generate the divergence matrix
void generateQT(int *QTRows, int *QTCols, real *QTVals, int nx, int ny);

// update the QT matrix for the method by Taira & Colonius
__global__ \
void updateQT(int *QTRows, int *QTCols, real *QTVals,
              int *ERows,  int *ECols,  real *EVals,
              int nx, int ny, real *x, real *y, real *dx,
              int totalPoints, real *xB, real *yB, int *I, int *J);

// update the QT matrix for the method by Taira & Colonius on the host
void updateQTHost(int *QTRows, int *QTCols, real *QTVals,
              int *ERows,  int *ECols,  real *EVals,
              int nx, int ny, real *x, real *y, real *dx,
              int totalPoints, real *xB, real *yB, int *I, int *J);

} // end of namespace kernels
