#pragma once

#include <types.h>
#include <helpers.h>

namespace kernels
{

__global__ \
void updateQ(int *QRows, int *QCols, real *QVals, int QSize, int *tags);

__global__ \
void updateQT(int *QTRows, int *QTCols, real *QTVals, int QTSize, int *tags, real *coeffs);

void generateQT(int *QTRows, int *QTCols, real *QTVals, int nx, int ny);

__global__ \
void updateQT(int *QTRows, int *QTCols, real *QTVals,
              int *ERows,  int *ECols,  real *EVals,
              int nx, int ny, real *x, real *y, real *dx,
              int totalPoints, real *xB, real *yB, int *I, int *J);

void updateQTHost(int *QTRows, int *QTCols, real *QTVals,
              int *ERows,  int *ECols,  real *EVals,
              int nx, int ny, real *x, real *y, real *dx,
              int totalPoints, real *xB, real *yB, int *I, int *J);

}
