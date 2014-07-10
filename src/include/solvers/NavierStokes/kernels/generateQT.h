/***************************************************************************//**
* \file generateQT.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the kernels required to generate matrix QT
*/

#pragma once

#include <types.h>
#include <helpers.h>

/********************//**
* \namespace kernels
* \brief Contain all custom-written CUDA kernels
*/
namespace kernels
{

/**
* \brief To be documented
*/
__global__ \
void updateQ(int *QRows, int *QCols, real *QVals, int QSize, int *tags);

/**
* \brief To be documented
*/
__global__ \
void updateQ(int *QRows, int *QCols, real *QVals, int QSize, int *tagsX, int *tagsY);

/**
* \brief To be documented
*/
void generateQT(int *QTRows, int *QTCols, real *QTVals, int nx, int ny);

/**
* \brief To be documented
*/
__global__ \
void updateQT(int *QTRows, int *QTCols, real *QTVals,
              int *ERows,  int *ECols,  real *EVals,
              int nx, int ny, real *x, real *y, real *dx,
              int totalPoints, real *xB, real *yB, int *I, int *J);

/**
* \brief to be documented
*/
void updateQTHost(int *QTRows, int *QTCols, real *QTVals,
              int *ERows,  int *ECols,  real *EVals,
              int nx, int ny, real *x, real *y, real *dx,
              int totalPoints, real *xB, real *yB, int *I, int *J);

} // end of namespace kernels
